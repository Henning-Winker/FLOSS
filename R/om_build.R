# Operating model: equilibrium calibration to MSY and stochastic hindcast
# driven by fleet effort.

#' Life-history FLPar for FLife
#' @param sp a [species()].
#' @return `FLPar` from [FLife::lhPar()].
#' @export
species_lhpar <- function(sp) {
  lh <- sp$lh
  FLife::lhPar(FLPar(linf = lh$linf, k = lh$k, t0 = lh$t0, a = lh$a,
                     b = lh$b, a50 = lh$tm, s = lh$s, m1 = lh$M,
                     tmax = lh$tmax))
}

#' Selectivity-at-length and -at-age by fleet
#'
#' Wraps [FLRef::sel_la()] for every fleet of the species, over ages
#' `0:tmax`.
#' @param sp a [species()].
#' @param lhpar optional `FLPar` from [species_lhpar()].
#' @return named list (one per fleet) of `sel_la()` results.
#' @export
species_sel <- function(sp, lhpar = species_lhpar(sp)) {
  lapply(sp$sel, function(z) {
    FLRef::sel_la(lhpar, amin = 0, amax = sp$lh$tmax, type = z$type,
                  L50 = z$L50, L95 = z$L95, Lpeak = z$Lpeak, sd = z$sd,
                  sd_left = z$sd_left, sd_right = z$sd_right)
  })
}

#' Equilibrium calibrated to a target MSY
#'
#' Builds the FLife equilibrium (`lhEql`, ages `0:tmax`, plus group `tmax`,
#' constant M), applies the joint selectivity-at-age (fleet selectivities
#' weighted by relative catchability, [FLRef::combine_sel_age()]) and
#' rescales virgin biomass so that equilibrium MSY equals `sp$msy` (two
#' passes of [FLRef::updsr()]).
#'
#' @param sp a [species()].
#' @param sels optional result of [species_sel()].
#' @return an `FLBRP` with reference points.
#' @export
calibrate_msy <- function(sp, sels = NULL) {
  lhpar <- species_lhpar(sp)
  if (is.null(sels)) sels <- species_sel(sp, lhpar)
  tmax <- sp$lh$tmax
  eql <- FLife::lhEql(lhpar, range = c(min = 0, max = tmax, minfbar = 1,
                                       maxfbar = tmax - 1, plusgroup = tmax))
  eql@m[] <- sp$lh$M
  jsel <- FLRef::combine_sel_age(sels, unlist(sp$relq))
  stk <- as(eql, "FLStock")
  units(stk) <- standardUnits(stk)
  sr <- as(as(eql, "predictModel"), "FLSR")
  stk@harvest[] <- c(jsel)
  out <- FLBRP::brp(FLBRP::FLBRP(stk, sr))
  for (i in 1:2) {
    v0 <- sp$msy / an(FLBRP::refpts(out)["msy", "yield"]) *
      an(FLBRP::refpts(out)["virgin", "ssb"])
    out <- FLBRP::brp(FLRef::updsr(out, s = sp$lh$s, v = v0))
  }
  out
}

#' Build a stochastic operating model for one species
#'
#' Apical fishing mortality of fleet \eqn{f} is
#' \eqn{F_{f,y} = q_f\,\eta_{f,y}\,E_{f,y}\,\delta_{f,y}} with
#' \eqn{q_f = relq_f \times F_{MSY}}, nominal effort \eqn{E} and efficiency
#' \eqn{\eta} from [effort_sim()] (shared by all species), and optional
#' species-specific catchability noise \eqn{\delta} (`sigmaQ`). Fleet F-at-age
#' is apical F times the fleet's selectivity-at-age; total F-at-age is the sum
#' over fleets. The stock is projected with [ffwd()] under Fbar
#' control and lognormal (optionally AR1) recruitment deviations. The first
#' two years are dropped as burn-in, as in the reef example.
#'
#' Reference points are recomputed for the realised exploitation pattern
#' (median F-at-age over iterations, last years) with
#' [FLRef::computeFbrp()] and stored in the returned `FLStockR`.
#'
#' @param sp a [species()].
#' @param fishery a [fishery()].
#' @param nsim number of iterations.
#' @param effort optional result of [effort_sim()] for this fishery (pass it
#'   when building several species so they share one effort realisation; by
#'   default it is regenerated from the fishery seed, which gives the same).
#' @param sigmaQ lognormal SD of species-specific catchability noise.
#' @param ar1 use the species' `rho` for AR1 recruitment deviations.
#' @param seed base seed for recruitment (and `sigmaQ`) deviations; combined
#'   with the species id so species differ but are reproducible.
#' @param treatment optional [treatment()] (changes to effort, recruitment,
#'   selectivity, creep or M); random numbers are unchanged, so a treated and
#'   an untreated operating model differ only by the treatment.
#' @return an object of class `floss_om`: a list with `stock` (`FLStockR`),
#'   `sr` (`FLSR`), `brp` (calibrated equilibrium), `sel` (fleet
#'   selectivities), `fapic` and `f_fleet` (`FLQuants` of apical F and
#'   F-at-age by fleet), `effort`, `efficiency`, `species`, `fishery` and
#'   `settings`.
#' @export
build_om <- function(sp, fishery, nsim = 100, effort = NULL, sigmaQ = 0,
                     ar1 = FALSE, seed = fishery$seed, treatment = NULL) {
  sp <- validate_species(sp, fishery)
  yrs <- fishery_years(fishery)
  ny <- length(yrs)
  if (is.null(effort)) effort <- effort_sim(fishery, nsim)
  if (dims(effort$effort[[1]])$iter != nsim)
    stop("'effort' has ", dims(effort$effort[[1]])$iter, " iterations, nsim = ",
         nsim)
  effort <- apply_effort_treatment(effort, treatment, yrs)

  lhpar <- species_lhpar(sp)
  sels <- species_sel(sp, lhpar)[fishery$fleets]
  eql <- calibrate_msy(sp, sels)
  fmsy <- an(FLBRP::refpts(eql)["msy", "harvest"])

  # starting stock: near-unfished, joint selectivity, all years
  start <- eql
  fbar(start) <- FLQuant(0.05 * fmsy, dimnames = list(year = yrs))
  stk <- as(start, "FLStock")
  units(stk) <- standardUnits(stk)
  sr <- as(as(eql, "predictModel"), "FLSR")
  stk <- propagate(stk, nsim)
  ages <- dimnames(stk)$age
  if (!is.null(treatment$M)) {
    yy <- ac(yrs[yrs >= treatment$M$from])
    m(stk)[, yy] <- m(stk)[, yy] * treatment$M$mult
  }

  sp_seed <- seed_for(seed, sp$id)
  res <- with_seed(sp_seed, {
    fapic <- lapply(fishery$fleets, function(f) {
      q <- sp$relq[[f]] * fmsy
      dq <- if (sigmaQ > 0) stats::rlnorm(ny * nsim, -0.5 * sigmaQ^2, sigmaQ)
            else 1
      out <- effort$effort[[f]] * effort$efficiency[[f]] * q
      out[] <- c(out) * dq
      units(out) <- "f"
      out
    })
    rho <- if (ar1) sp$lh$rho else 0
    devs <- rlnormar1(n = nsim, meanlog = 0, sdlog = sp$lh$sigR,
                             rho = rho, years = yrs[-1], bias.correct = TRUE)
    devs <- devs * rec_multiplier(treatment, yrs[-1])
    list(fapic = fapic, devs = devs)
  })
  fapic <- FLQuants(stats::setNames(res$fapic, fishery$fleets))

  # fleet F-at-age = apical F x selectivity-at-age
  f_fleet <- FLQuants(lapply(fishery$fleets, function(f) {
    s_a <- c(sels[[f]]$sel_a)
    sa_y <- matrix(s_a, length(ages), length(yrs))
    for (sc in Filter(function(z) identical(z$fleet, f), treatment$selectivity)) {
      new <- FLRef::sel_la(lhpar, amin = 0, amax = sp$lh$tmax, type = sc$sel$type,
                           L50 = sc$sel$L50, L95 = sc$sel$L95, Lpeak = sc$sel$Lpeak,
                           sd = sc$sel$sd, sd_left = sc$sel$sd_left,
                           sd_right = sc$sel$sd_right)
      sa_y[, yrs >= sc$from] <- c(new$sel_a)
    }
    flq(rep(c(fapic[[f]]), each = length(ages)) * c(sa_y), age = ages,
        years = yrs, nsim = nsim, units = "f")
  }))
  names(f_fleet) <- fishery$fleets
  ftot <- Reduce(`+`, f_fleet)
  harvest(stk)[, -1] <- ftot[, -1]
  units(harvest(stk)) <- "f"
  fb <- fbar(stk)

  om <- ffwd(stk, sr = sr, deviances = res$devs,
                      control = FLasher::fwdControl(
                        year = yrs[-1], quant = "fbar",
                        value = fb[, ac(yrs[-1])]))
  om <- window(om, start = yrs[3])
  landings(om) <- computeLandings(om)
  discards(om) <- computeDiscards(om)
  catch(om) <- computeCatch(om)

  keep <- ac(yrs[3]:yrs[ny])
  f_fleet <- FLQuants(lapply(f_fleet, function(x) x[, keep]))
  fapic <- FLQuants(lapply(fapic, function(x) x[, keep]))

  # reference points for the realised exploitation pattern
  ref_stk <- iter(om, 1)
  harvest(ref_stk) <- iterMedians(harvest(om))
  eq <- FLRef::computeFbrp(ref_stk, sr, proxy = "msy", blim = 0.25,
                           type = "btgt", verbose = FALSE)
  pts <- FLRef::Fbrp(eq)
  bmsy <- an(pts[2])
  om <- FLStockR(om)
  om@refpts <- FLPar(Fmsy = an(pts[1]), Bmsy = bmsy, Bthr = 0.5 * bmsy,
                     Blim = 0.25 * bmsy,
                     MSY = an(FLBRP::refpts(eq)["msy", "yield"]),
                     B0 = an(FLBRP::refpts(eq)["virgin", "ssb"]))
  name(om) <- sp$species
  desc(om) <- paste("FLOSS OM:", sp$species, "|", fishery$name)

  structure(list(
    stock = om, sr = sr, brp = eql, brp_realised = eq, sel = sels,
    fapic = fapic, f_fleet = f_fleet,
    effort = FLQuants(lapply(effort$effort, function(x) x[, keep])),
    efficiency = FLQuants(lapply(effort$efficiency, function(x) x[, keep])),
    species = sp, fishery = fishery, treatment = treatment,
    settings = list(nsim = nsim, sigmaQ = sigmaQ, ar1 = ar1, seed = seed,
                    species_seed = sp_seed, fishery_hash = spec_hash(fishery),
                    floss_version = as.character(utils::packageVersion("FLOSS")),
                    built = Sys.time())),
    class = "floss_om")
}

#' Build operating models for several species sharing one fishery
#'
#' All species see the same effort realisation; recruitment deviations
#' differ by species. Pass `oms` to reuse operating models already built for
#' the same fishery (matched by species id and fishery hash), so adding a
#' species does not rebuild the others.
#'
#' @param species_list a list of [species()] (or species ids from the
#'   library).
#' @inheritParams build_om
#' @param oms optional named list of existing `floss_om` objects to reuse.
#' @param ... passed to [build_om()].
#' @return a named list of class `floss_oms`.
#' @export
build_oms <- function(species_list, fishery, nsim = 100, oms = NULL, ...) {
  species_list <- lapply(species_list, function(s)
    if (is.character(s)) load_species(s) else s)
  effort <- effort_sim(fishery, nsim)
  h <- spec_hash(fishery)
  out <- lapply(species_list, function(s) {
    old <- oms[[s$id]]
    if (!is.null(old) && identical(old$settings$fishery_hash, h) &&
        identical(old$settings$nsim, nsim) &&
        identical(spec_hash(old$species), spec_hash(s))) return(old)
    build_om(s, fishery, nsim = nsim, effort = effort, ...)
  })
  names(out) <- vapply(species_list, `[[`, "", "id")
  structure(out, class = "floss_oms")
}

#' True catch by fleet
#'
#' Partial Baranov catch of each fleet from fleet F-at-age, total Z and
#' numbers-at-age of the operating model:
#' \eqn{C_{f,a,y} = N_{a,y}\,F_{f,a,y}/Z_{a,y}\,(1 - e^{-Z_{a,y}})}.
#'
#' @param om a `floss_om`.
#' @param type `"wt"` (catch weight, t) or `"n"` (catch-at-age, numbers).
#' @return `FLQuants` by fleet.
#' @export
catch_fleet <- function(om, type = c("wt", "n")) {
  type <- match.arg(type)
  stk <- om$stock
  z <- harvest(stk) + m(stk)
  surv <- stock.n(stk) * (1 - exp(-z)) / z
  out <- lapply(om$f_fleet, function(f) {
    cn <- surv * f
    if (type == "n") cn else quantSums(cn * catch.wt(stk))
  })
  FLQuants(out)
}

#' @exportS3Method base::print
print.floss_om <- function(x, ...) {
  s <- x$stock
  cat("FLOSS operating model:", x$species$species, "\n")
  cat("  fishery:", x$fishery$name, "(", toString(x$fishery$fleets), ")\n")
  cat("  years:", dims(s)$minyear, "-", dims(s)$maxyear, " iterations:",
      dims(s)$iter, "\n")
  rp <- x$stock@refpts
  cat("  refpts:", paste(rownames(rp), signif(c(rp), 3), sep = "=",
                         collapse = ", "), "\n")
  last <- ac(dims(s)$maxyear)
  cat("  final year (median): B/Bmsy =",
      signif(stats::median(c(ssb(s)[, last])) / c(rp["Bmsy"]), 3),
      " F/Fmsy =", signif(stats::median(c(fbar(s)[, last])) / c(rp["Fmsy"]), 3),
      "\n")
  invisible(x)
}

spec_hash <- function(x) {
  f <- tempfile()
  on.exit(unlink(f))
  saveRDS(unclass_deep(x), f, version = 2, compress = FALSE)
  unname(tools::md5sum(f))
}
