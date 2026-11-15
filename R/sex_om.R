# Sex-structured operating model: females and males in the unit dimension of
# an FLStock (units "F", "M"), projected with FLasher.
#
# FLasher with units: recruitment follows the summed spawning output of all
# units (male maturity is set to 0, so this is female SSB) and recruits are
# split equally between the units (sex ratio 0.5 at recruitment). Fishing is
# controlled by 'effort': F-at-age = effort x harvest pattern, with the
# pattern scaled to a maximum of 1 over both sexes, so effort is the apical F
# at length shared by the sexes (the F of length-based methods such as FLicc).
# FLBRP does not handle units, so reference points are computed here from the
# sex-structured equilibrium (refpts_sex()) and can be checked against
# deterministic FLasher projections (refpts_fwd_sex()).

#' Life history by sex
#'
#' Female and male life-history parameters for a sex-structured operating
#' model. Male parameters not given are copied from the females.
#'
#' @param female named list: `linf`, `k`, `t0` (von Bertalanffy, cm and
#'   years), `a`, `b` (W = aL^b), `L50`, optional `L95` (maturity at length;
#'   default `L50 + 0.2 linf`, as FLicc), `M`, optional `cvl` (CV of length
#'   at age, default 0.1).
#' @param male named list of the male values that differ.
#' @param ... male values as named arguments, alternative to `male`.
#' @return a list of class `floss_lh_sex` with elements `F` and `M`.
#' @export
#' @examples
#' lh <- lh_sex(list(linf = 60, k = 0.15, t0 = -0.3, a = 1e-5, b = 3, L50 = 35, M = 0.3),
#'              linf = 45, k = 0.25, L50 = 28)
lh_sex <- function(female, male = NULL, ...) {
  dots <- list(...)
  if (!is.null(male) && length(dots)) stop("Give male values in 'male' or as arguments, not both.")
  male <- if (is.null(male)) dots else as.list(male)
  female <- as.list(female)
  need <- c("linf", "k", "t0", "a", "b", "L50", "M")
  miss <- setdiff(need, names(female))
  if (length(miss)) stop("missing female parameter(s): ", toString(miss))
  fill <- function(x) {
    if (is.null(x$cvl)) x$cvl <- 0.1
    if (is.null(x$L95)) x$L95 <- x$L50 + 0.2 * x$linf
    lapply(x, as.numeric)
  }
  fem <- fill(female)
  mal <- female; mal[names(male)] <- male
  if (is.null(male$L95)) mal$L95 <- NULL
  mal <- fill(mal)
  structure(list(F = fem, M = mal), class = "floss_lh_sex")
}

# selectivity at length for a floss_sel (also accepts FLicc's "dsnormal")
sel_len_floss <- function(sel, L) {
  type <- if (sel$type == "dsnormal") "dnormal" else sel$type
  switch(type,
    logistic = 1 / (1 + exp(-log(19) * (L - sel$L50) / (sel$L95 - sel$L50))),
    normal = exp(-(L - sel$Lpeak)^2 / (2 * sel$sd^2)),
    dnormal = ifelse(L < sel$Lpeak, exp(-(L - sel$Lpeak)^2 / (2 * sel$sd_left^2)),
                     exp(-(L - sel$Lpeak)^2 / (2 * sel$sd_right^2))),
    stop("unknown selectivity type ", sel$type))
}

#' Age-length keys by sex
#'
#' Probability of each length bin given age, for one sex: normal length at
#' age around the von Bertalanffy mean at `age + timing` with SD
#' `cvl x mean`.
#'
#' @param lh one sex of [lh_sex()].
#' @param ages ages.
#' @param lens lower bounds of the length bins (equal width).
#' @param timing fraction of the year (0 start, 0.5 mid-year).
#' @return matrix ages x length bins (rows sum to 1).
#' @export
alk_floss <- function(lh, ages, lens, timing = 0) {
  w <- if (length(lens) > 1) lens[2] - lens[1] else 1
  mu <- pmax(lh$linf * (1 - exp(-lh$k * (ages + timing - lh$t0))), 1e-3)
  sd <- pmax(lh$cvl * mu, 1e-3)
  ub <- c(lens[-1], Inf); lb <- c(-Inf, lens[-1])
  p <- t(sapply(seq_along(ages), function(i)
    stats::pnorm(ub, mu[i], sd[i]) - stats::pnorm(lb, mu[i], sd[i])))
  dimnames(p) <- list(age = ages, len = lens)
  p / rowSums(p)
}

#' Sex-structured FLStock
#'
#' Builds an `FLStock` with units `"F"` and `"M"`: natural mortality,
#' maturity (female; males 0, so `ssb()` and the stock-recruitment
#' relationship use female spawning biomass), stock weights at the start of
#' the year and catch weights at mid-year from the length-at-age
#' distributions, and the harvest pattern of the fleets (selectivity at
#' length integrated over length at mid-year, weighted by relative
#' catchability and scaled to a maximum of 1 over both sexes).
#'
#' @param lh [lh_sex()].
#' @param sel named list of [selectivity()] objects, one per fleet.
#' @param relq relative catchability by fleet (default equal).
#' @param tmax maximum age (plus group).
#' @param years years.
#' @param lens lower bounds of the length bins (cm).
#' @return an `FLStock` with attributes `sel_fleet` (`FLQuants` of fleet
#'   selectivity at age by sex, on the scale of the joint pattern), `lh`,
#'   `lens` and `sel`.
#' @export
stock_sex <- function(lh, sel, relq = NULL, tmax = 20, years = 2000:2025,
                      lens = seq(0, ceiling(1.3 * max(lh$F$linf, lh$M$linf)), by = 2)) {
  ages <- 0:tmax
  sel <- lapply(sel, function(z) if (inherits(z, "floss_sel")) z else do.call(selectivity, z))
  if (is.null(relq)) relq <- stats::setNames(rep(1, length(sel)), names(sel))
  relq <- relq[names(sel)]
  mid <- lens + diff(lens[1:2]) / 2
  dmn <- list(age = ages, year = years, unit = c("F", "M"))
  q <- function(v) FLQuant(v, dimnames = dmn)
  m <- mat <- swt <- cwt <- q(0)
  sfl <- lapply(sel, function(z) q(0))
  for (s in c("F", "M")) {
    l <- lh[[s]]
    a0 <- alk_floss(l, ages, lens, 0); a5 <- alk_floss(l, ages, lens, 0.5)
    matL <- 1 / (1 + exp(-log(19) * (mid - l$L50) / (l$L95 - l$L50)))
    m[, , s] <- l$M
    mat[, , s] <- if (s == "F") c(a0 %*% matL) else 0
    swt[, , s] <- c(a0 %*% (l$a * mid^l$b))
    cwt[, , s] <- c(a5 %*% (l$a * mid^l$b))
    for (f in names(sel)) sfl[[f]][, , s] <- relq[f] * c(a5 %*% sel_len_floss(sel[[f]], mid))
  }
  joint <- Reduce(`+`, sfl)
  mx <- max(joint)
  sfl <- lapply(sfl, function(x) x / mx)
  stk <- FLStock(stock.n = q(1), m = m, mat = mat, stock.wt = swt,
                 catch.wt = cwt, landings.wt = cwt, discards.wt = cwt,
                 harvest = joint / mx)
  units(harvest(stk)) <- "f"
  harvest.spwn(stk)[] <- 0; m.spwn(stk)[] <- 0
  discards.n(stk)[] <- 0; landings.n(stk)[] <- 1; catch.n(stk)[] <- 1
  range(stk, c("minfbar", "maxfbar")) <- range(ages[ages > 0 & ages < tmax])
  attr(stk, "sel_fleet") <- FLQuants(sfl)
  attr(stk, "lh") <- lh; attr(stk, "lens") <- lens; attr(stk, "sel") <- sel
  stk
}

# year-averaged biology and harvest pattern by sex (last nyears)
eq_inputs_sex <- function(stk, nyears = 3) {
  yy <- utils::tail(dimnames(stk)$year, nyears)
  g <- function(x) apply(x@.Data[, yy, , 1, 1, 1, drop = FALSE], c(1, 3), mean)
  h <- g(harvest(stk)); h <- h / max(h)
  list(m = g(m(stk)), mat = g(mat(stk)), swt = g(stock.wt(stk)),
       cwt = g(catch.wt(stk)), sel = h)
}

# numbers, catch per recruit of one sex at apical F (plus group in last age)
pr_sex1 <- function(F, m, sel, cwt) {
  z <- m + F * sel; na <- length(z)
  n <- numeric(na); n[1] <- 1
  for (a in 2:na) n[a] <- n[a - 1] * exp(-z[a - 1])
  n[na] <- n[na] / (1 - exp(-z[na]))
  cn <- n * F * sel / z * (1 - exp(-z))
  list(n = n, cn = cn, yw = sum(cn * cwt))
}

#' Equilibrium of a sex-structured stock
#'
#' Per-recruit and equilibrium quantities by apical F for an `FLStock` with
#' units `"F"`, `"M"` ([stock_sex()]) and a Beverton-Holt
#' stock-recruitment relationship ([sr_sex()]). Recruits are split equally
#' between the sexes (as FLasher); spawning biomass is female.
#'
#' @param stk sex-structured `FLStock`.
#' @param sr `predictModel` or `FLSR` (Beverton-Holt, `a`, `b`), or `NULL`
#'   for per-recruit quantities only.
#' @param F apical F values.
#' @param nyears years averaged for biology and harvest pattern.
#' @return data.frame with `F`, `spr` (female SSB per recruit relative to
#'   unfished), `sbpr`, `ypr`, `rec`, `ssb`, `yield` and `yield_F`,
#'   `yield_M` (yield by sex).
#' @export
eq_sex <- function(stk, sr = NULL, F = seq(0, 2, by = 0.01), nyears = 3) {
  x <- eq_inputs_sex(stk, nyears)
  sb <- function(f) 0.5 * sum(pr_sex1(f, x$m[, "F"], x$sel[, "F"], x$cwt[, "F"])$n *
                              x$mat[, "F"] * x$swt[, "F"])
  sb0 <- sb(0)
  out <- do.call(rbind, lapply(F, function(f) {
    pF <- pr_sex1(f, x$m[, "F"], x$sel[, "F"], x$cwt[, "F"])
    pM <- pr_sex1(f, x$m[, "M"], x$sel[, "M"], x$cwt[, "M"])
    sbpr <- 0.5 * sum(pF$n * x$mat[, "F"] * x$swt[, "F"])
    data.frame(F = f, sbpr = sbpr, spr = sbpr / sb0, ypr = 0.5 * (pF$yw + pM$yw),
               ypr_F = 0.5 * pF$yw, ypr_M = 0.5 * pM$yw)
  }))
  if (!is.null(sr)) {
    p <- params(sr)
    a <- an(p["a"]); b <- an(p["b"])
    out$rec <- pmax(a - b / out$sbpr, 0)
    out$ssb <- out$rec * out$sbpr
    out$yield <- out$rec * out$ypr
    out$yield_F <- out$rec * out$ypr_F; out$yield_M <- out$rec * out$ypr_M
  }
  out
}

#' Beverton-Holt stock-recruitment for a sex-structured stock
#'
#' Steepness parameterisation on female spawning biomass, with unfished
#' recruitment `R0` (both sexes) and unfished female SSB per recruit from
#' the stock's biology.
#'
#' @param stk sex-structured `FLStock`.
#' @param s steepness.
#' @param R0 unfished recruitment (both sexes); or give `B0`.
#' @param B0 unfished female spawning biomass (alternative to `R0`).
#' @return a `predictModel` (`bevholt`, `a`, `b`) for FLasher, with
#'   attribute `pars` (`s`, `R0`, `B0`).
#' @export
sr_sex <- function(stk, s = 0.7, R0 = 1000, B0 = NULL) {
  sbpr0 <- eq_sex(stk, F = 0)$sbpr
  if (!is.null(B0)) R0 <- B0 / sbpr0
  v <- R0 * sbpr0
  a <- 4 * s * R0 / (5 * s - 1); b <- v * (1 - s) / (5 * s - 1)
  # only a and b in the model (FLasher reads all parameters of the model);
  # steepness, R0 and B0 kept as an attribute
  out <- predictModel(model = bevholt()$model, params = FLPar(a = a, b = b))
  attr(out, "pars") <- c(s = s, R0 = R0, B0 = v)
  out
}

#' Reference points of a sex-structured stock
#'
#' MSY reference points (yield of both sexes maximised over apical F, female
#' spawning biomass), unfished levels, Fcrash and F at given female SPR,
#' from [eq_sex()].
#'
#' @inheritParams eq_sex
#' @param spr female SPR levels for Fspr (e.g. 0.4).
#' @return `FLPar` with `Fmsy`, `Bmsy`, `MSY`, `B0`, `R0`, `Fcrash`,
#'   `SPRmsy`, `Bmsy_B0` and `Fspr<x>`.
#' @export
refpts_sex <- function(stk, sr, spr = c(0.3, 0.4), nyears = 3) {
  yf <- function(f) eq_sex(stk, sr, F = f, nyears)$yield
  sb <- function(f) eq_sex(stk, sr, F = f, nyears)$ssb
  e0 <- eq_sex(stk, sr, F = 0, nyears)
  fhi <- 1; while (sb(fhi) > 0 && fhi < 50) fhi <- fhi * 2
  fcrash <- stats::uniroot(function(f) sb(f) - 1e-6 * e0$ssb, c(0, fhi))$root
  fmsy <- stats::optimize(yf, c(0, fcrash), maximum = TRUE, tol = 1e-8)$maximum
  em <- eq_sex(stk, sr, F = fmsy, nyears)
  sprf <- function(f) eq_sex(stk, NULL, F = f, nyears)$spr
  fs <- sapply(spr, function(p) stats::uniroot(function(f) sprf(f) - p, c(0, 50))$root)
  out <- c(Fmsy = fmsy, Bmsy = em$ssb, MSY = em$yield, B0 = e0$ssb, R0 = e0$rec,
           Fcrash = fcrash, SPRmsy = em$spr, Bmsy_B0 = em$ssb / e0$ssb)
  out <- c(out, stats::setNames(fs, paste0("Fspr", round(100 * spr))))
  do.call(FLPar, as.list(out))
}

#' MSY by deterministic projection (check of refpts_sex())
#'
#' Projects the sex-structured stock with FLasher at constant apical F
#' (`effort`) for `nyears` without recruitment deviations, and maximises the
#' final-year catch over F (golden-section search). Slow; for checking the
#' equilibrium reference points.
#'
#' @inheritParams eq_sex
#' @param nyears length of the projection.
#' @param interval search interval for F.
#' @return `FLPar` with `Fmsy`, `Bmsy`, `MSY`.
#' @export
refpts_fwd_sex <- function(stk, sr, nyears = 100, interval = c(0.01, 2)) {
  # start from a positive F: FLasher limits effort to 100 times the
  # effort of the year before the projection
  st <- init_stock_sex(stk, sr, F = 0.05, years = seq(2000, length.out = nyears + 1))
  proj <- function(f) {
    o <- FLasher::fwd(st, sr = sr, control = FLasher::fwdControl(
      year = 2000 + seq_len(nyears), quant = "effort", value = f))
    yr <- ac(2000 + nyears)
    c(catch = sum(catch.n(o)[, yr] * catch.wt(o)[, yr]), ssb = sum(ssb(o)[, yr]))
  }
  opt <- stats::optimize(function(f) proj(f)["catch"], interval, maximum = TRUE, tol = 1e-5)
  r <- proj(opt$maximum)
  FLPar(Fmsy = opt$maximum, Bmsy = an(r["ssb"]), MSY = an(r["catch"]))
}

# stock over 'years' in equilibrium at apical F in the first year
init_stock_sex <- function(stk, sr, F, years) {
  st <- expand_years_sex(stk, years)
  x <- eq_inputs_sex(stk)
  e <- eq_sex(stk, sr, F = F)
  for (s in c("F", "M")) {
    p <- pr_sex1(F, x$m[, s], x$sel[, s], x$cwt[, s])
    stock.n(st)[, 1, s] <- 0.5 * e$rec * p$n
  }
  harvest(st)[, 1] <- F * harvest(st)[, 1] / max(harvest(st)[, 1])
  st
}

# FLStock over new years with the last-year biology and harvest pattern
expand_years_sex <- function(stk, years) {
  ly <- utils::tail(dimnames(stk)$year, 1)
  one <- stk[, ly]
  dmn <- dimnames(stock.n(stk)); dmn$year <- ac(years)
  rep_y <- function(x) FLQuant(rep(c(x), length(years)), dimnames = dmn, units = units(x))
  rep_y2 <- function(x) {
    a <- array(c(x), dim = c(dim(x)[1], 1, dim(x)[3], 1, 1, dim(x)[6]))
    out <- FLQuant(NA, dimnames = c(dmn[1:5], list(iter = dimnames(x)$iter)), units = units(x))
    for (i in seq_along(years)) out[, i] <- a
    out
  }
  st <- FLStock(stock.n = rep_y2(stock.n(one)), m = rep_y2(m(one)), mat = rep_y2(mat(one)),
                stock.wt = rep_y2(stock.wt(one)), catch.wt = rep_y2(catch.wt(one)),
                landings.wt = rep_y2(landings.wt(one)), discards.wt = rep_y2(discards.wt(one)),
                harvest = rep_y2(harvest(one)))
  harvest.spwn(st)[] <- 0; m.spwn(st)[] <- 0
  discards.n(st)[] <- 0; landings.n(st)[] <- 1; catch.n(st)[] <- 1
  units(harvest(st)) <- "f"
  range(st, c("minfbar", "maxfbar")) <- range(stk)[c("minfbar", "maxfbar")]
  for (a in c("sel_fleet", "lh", "lens", "sel")) attr(st, a) <- attr(stk, a)
  st
}

#' Build a sex-structured operating model
#'
#' Starts from the equilibrium at the first year's apical F and projects the
#' stock with FLasher under apical F (`effort`) with lognormal (optionally
#' AR1) recruitment deviations. Recruits are split equally between the
#' sexes; spawning biomass is female.
#'
#' @param stk sex-structured `FLStock` from [stock_sex()].
#' @param sr stock-recruitment from [sr_sex()].
#' @param fapic apical F by year (numeric vector over `years`, or `FLQuant`).
#' @param years years (default the stock's years).
#' @param nsim iterations.
#' @param sigR recruitment SD (log scale).
#' @param rho recruitment autocorrelation.
#' @param seed random seed.
#' @return a list of class `floss_om_sex`: `stock` (`FLStock` with units F,
#'   M), `sr`, `refpts` ([refpts_sex()]), `fapic` and `truth` (by year and
#'   iteration: apical F, F/Fmsy, female SSB, SSB/B0, SSB/Bmsy, female SPR
#'   at the year's F, catch).
#' @export
om_sex <- function(stk, sr, fapic, years = as.numeric(dimnames(stk)$year),
                   nsim = 1, sigR = 0.5, rho = 0, seed = 1) {
  fapic <- as.numeric(fapic)
  if (length(fapic) != length(years)) stop("fapic needs one value per year.")
  st <- init_stock_sex(stk, sr, F = max(fapic[1], 1e-3), years = years)
  st <- propagate(st, nsim)
  set.seed(seed)
  ny <- length(years) - 1
  e <- matrix(stats::rnorm(ny * nsim, 0, sigR), ny, nsim)
  if (rho > 0) for (i in 2:ny) e[i, ] <- rho * e[i - 1, ] + sqrt(1 - rho^2) * e[i, ]
  devs <- FLQuant(c(exp(e - sigR^2 / 2)), dimnames = list(age = 0, year = years[-1], iter = seq_len(nsim)))
  om <- FLasher::fwd(st, sr = sr, deviances = devs,
                     control = FLasher::fwdControl(year = years[-1], quant = "effort", value = fapic[-1]))
  landings.n(om) <- catch.n(om); landings(om) <- computeLandings(om)
  discards(om) <- computeDiscards(om); catch(om) <- computeCatch(om)
  rp <- refpts_sex(stk, sr)
  out <- structure(list(stock = om, sr = sr, refpts = rp, fapic = fapic, years = years),
                   class = "floss_om_sex")
  out$truth <- truth_sex(out)
  out
}

#' True status of a sex-structured operating model
#'
#' @param om [om_sex()] result.
#' @return data.frame by year and iteration: `F` (apical), `F_Fmsy`, `ssb`
#'   (female), `ssb_B0`, `ssb_Bmsy`, `spr` (female SPR at the year's F,
#'   equilibrium per recruit), `catch`, `rec`.
#' @export
truth_sex <- function(om) {
  st <- om$stock; rp <- om$refpts
  sprF <- eq_sex(st, NULL, F = om$fapic, nyears = 1)$spr
  ssbq <- unitSums(ssb(st))
  it <- dimnames(st)$iter
  do.call(rbind, lapply(seq_along(it), function(i) {
    s <- c(ssbq[, , , , , i])
    data.frame(year = om$years, iter = it[i], F = om$fapic, F_Fmsy = om$fapic / an(rp["Fmsy"]),
               ssb = s, ssb_B0 = s / an(rp["B0"]), ssb_Bmsy = s / an(rp["Bmsy"]), spr = sprF,
               catch = c(unitSums(catch(st))[, , , , , i]),
               rec = c(unitSums(stock.n(st)[1])[, , , , , i]))
  }))
}

#' Length frequencies by sex from a sex-structured operating model
#'
#' Expected numbers at length by sex of a fleet's catch (fleet catch at age
#' from the share of the fleet in F at age, through mid-year age-length
#' keys) or of a survey (numbers at age at `timing` times the survey's
#' selectivity at length), and a multinomial sample of `n` fish per year.
#' Sexed samples are drawn jointly over length and sex, so the sex ratio at
#' length is kept; unsexed samples are the sum over sexes.
#'
#' @param om [om_sex()] result.
#' @param fleet name of a fleet of the stock (catch samples), or `NULL` for
#'   a survey.
#' @param survey_sel a [selectivity()] for a survey (when `fleet = NULL`).
#' @param n fish sampled per year (`Inf` for expected proportions).
#' @param sexed `TRUE` for samples by sex (units F, M), `FALSE` for unsexed
#'   (unit U).
#' @param years years sampled (default all).
#' @param timing fraction of the year of the survey.
#' @param iter iteration.
#' @return `FLQuant` length x year x unit, ready for FLicc.
#' @export
lfd_sex <- function(om, fleet = NULL, survey_sel = NULL, n = 400, sexed = FALSE,
                    years = om$years, timing = 0.5, iter = 1) {
  st <- om$stock; lh <- attr(st, "lh"); lens <- attr(st, "lens")
  ages <- as.numeric(dimnames(st)$age); mid <- lens + diff(lens[1:2]) / 2
  yy <- ac(years)
  ex <- lapply(c(F = "F", M = "M"), function(s) {
    if (!is.null(fleet)) {
      sf <- attr(st, "sel_fleet")[[fleet]][, 1, s, , , 1, drop = TRUE]
      f_fl <- outer(c(sf), om$fapic[match(years, om$years)])
      fz <- f_fl / (harvest(st)[, yy, s, , , iter, drop = TRUE] + m(st)[, yy, s, , , iter, drop = TRUE])
      z <- harvest(st)[, yy, s, , , iter, drop = TRUE] + m(st)[, yy, s, , , iter, drop = TRUE]
      ca <- stock.n(st)[, yy, s, , , iter, drop = TRUE] * fz * (1 - exp(-z))
      t(alk_floss(lh[[s]], ages, lens, 0.5)) %*% ca
    } else {
      z <- harvest(st)[, yy, s, , , iter, drop = TRUE] + m(st)[, yy, s, , , iter, drop = TRUE]
      na <- stock.n(st)[, yy, s, , , iter, drop = TRUE] * exp(-timing * z)
      (t(alk_floss(lh[[s]], ages, lens, timing)) %*% na) * sel_len_floss(survey_sel, mid)
    }
  })
  units <- if (sexed) c("F", "M") else "U"
  out <- FLQuant(0, dimnames = list(len = lens, year = yy, unit = units), units = "cm")
  for (j in seq_along(yy)) {
    pF <- ex$F[, j]; pM <- ex$M[, j]
    if (sexed) {
      p <- c(pF, pM); x <- if (is.finite(n)) stats::rmultinom(1, n, p / sum(p))[, 1] else p / sum(p)
      out[, j, "F"] <- x[seq_along(lens)]; out[, j, "M"] <- x[length(lens) + seq_along(lens)]
    } else {
      p <- pF + pM
      out[, j, "U"] <- if (is.finite(n)) stats::rmultinom(1, n, p / sum(p))[, 1] else p / sum(p)
    }
  }
  out
}

#' Life history by sex for FLicc
#'
#' @param lh [lh_sex()].
#' @return `FLPar` params x unit (F, M) with `linf`, `k`, `t0`, `M`, `L50`,
#'   `L95`, `a`, `b`, as used by `FLicc::stocklen_sex()`.
#' @export
lhpar_flicc_sex <- function(lh) {
  p <- c("linf", "k", "t0", "M", "L50", "L95", "a", "b")
  a <- array(c(unlist(lh$F[p]), unlist(lh$M[p])), dim = c(length(p), 2, 1),
             dimnames = list(params = p, unit = c("F", "M"), iter = 1))
  FLPar(a)
}
