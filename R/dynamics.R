# Views of the true dynamics of an operating model: equilibrium curves with
# the realised trajectory, stock-recruitment, and the effect of selectivity.

#' Equilibrium curves with the realised trajectory
#'
#' Yield, SSB and recruitment against F and SSB at equilibrium
#' ([FLRef::ploteq()]), with the stock's yearly values as points
#' (`obs = TRUE`). During a fishing-down phase catches lie above the
#' equilibrium yield curve: high catches are possible only while biomass is
#' being removed. The curves are for the exploitation pattern of the last
#' years; with a changing fleet mix, earlier years belonged to a different
#' equilibrium.
#'
#' @param om a `floss_om`.
#' @param iter iteration whose trajectory is shown.
#' @param obs show the realised trajectory.
#' @param rel plot relative to the MSY reference points.
#' @return a `ggplot`.
#' @export
plot_eq <- function(om, iter = 1, obs = TRUE, rel = FALSE) {
  eq <- FLRef::computeFbrp(FLCore::iter(om$stock, iter), om$sr, proxy = "msy",
                           blim = 0.25, type = "btgt", verbose = FALSE)
  FLRef::ploteq(eq, obs = obs, rel = rel) +
    ggplot2::ggtitle(om$species$species,
                     subtitle = if (obs) paste("Points: realised years, iteration", iter))
}

#' Stock-recruitment curve with simulated pairs
#'
#' The Beverton-Holt curve of the operating model with the SSB-recruitment
#' pairs of the selected iterations (recruitment at age 0 in the same year).
#' The scatter is the recruitment variability (sigmaR); it is why steepness
#' is hard to estimate from data.
#' @param om a `floss_om`.
#' @param iters iterations to show.
#' @return a `ggplot`.
#' @export
plot_srr <- function(om, iters = 1:3) {
  s <- FLCore::iter(om$stock, iters)
  pr <- params(om$sr)
  a <- c(pr["a"])
  b <- c(pr["b"])
  d <- data.frame(ssb = c(ssb(s)), rec = c(rec(s)),
                  iter = factor(rep(iters, each = dims(s)$year)))
  S <- seq(0, max(d$ssb, c(om$stock@refpts["B0"])) * 1.05, length.out = 200)
  ggplot2::ggplot(d, ggplot2::aes(.data$ssb, .data$rec)) +
    ggplot2::geom_line(data = data.frame(ssb = S, rec = a * S / (b + S)),
                       colour = "#0A5F78", linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$iter), alpha = 0.6, size = 1.6) +
    ggplot2::geom_vline(xintercept = c(om$stock@refpts["Bmsy"]), linetype = 2) +
    ggplot2::labs(x = "SSB", y = "Recruitment (age 0)", colour = "Iteration",
                  title = om$species$species,
                  subtitle = sprintf("Beverton-Holt, steepness %.2f, sigmaR %.2f; dashed: Bmsy",
                                     om$species$lh$s, om$species$lh$sigR)) +
    ggplot2::theme_bw()
}

#' Equilibrium curves under alternative selectivities
#'
#' The same stock (life history, recruitment, unfished biomass) fished with
#' different selectivities: equilibrium yield and SSB against F (Fbar).
#' Catching fish before they have grown lowers the maximum yield and the F
#' at which it is reached (growth overfishing); protecting small fish raises
#' both. Points mark MSY for each selectivity.
#'
#' @param sp a [species()].
#' @param sels named list of [selectivity()] objects.
#' @param fmax largest F shown, as a multiple of Fmsy of the species' own
#'   (joint) selectivity.
#' @return a `ggplot`.
#' @export
#' @examples
#' \donttest{
#' sp <- load_species("lutjanus_bohar")
#' plot_eq_selectivity(sp, list(
#'   small = selectivity("logistic", L50 = 25, L95 = 30),
#'   mature = selectivity("logistic", L50 = 35, L95 = 42),
#'   large = selectivity("logistic", L50 = 50, L95 = 58)))
#' }
plot_eq_selectivity <- function(sp, sels, fmax = 10) {
  base <- calibrate_msy(sp)                 # stock scaled to the target MSY
  fgrid <- seq(0, fmax * an(FLBRP::refpts(base)["msy", "harvest"]), length.out = 101)
  lhpar <- species_lhpar(sp)
  curves <- do.call(rbind, lapply(names(sels), function(n) {
    z <- sels[[n]]
    sa <- FLRef::sel_la(lhpar, amin = 0, amax = sp$lh$tmax, type = z$type,
                        L50 = z$L50, L95 = z$L95, Lpeak = z$Lpeak, sd = z$sd,
                        sd_left = z$sd_left, sd_right = z$sd_right)$sel_a
    b <- base
    landings.sel(b)[] <- c(sa)              # same stock, new selectivity
    discards.sel(b)[] <- 0
    fbar(b) <- FLQuant(fgrid)
    b <- FLBRP::brp(b)
    data.frame(sel = n, F = fgrid, yield = c(landings(b)), ssb = c(ssb(b)))
  }))
  curves$sel <- factor(curves$sel, levels = names(sels))
  msy <- do.call(rbind, lapply(split(curves, curves$sel), function(d) d[which.max(d$yield), ]))
  long <- rbind(data.frame(curves[, c("sel", "F")], what = "Yield", value = curves$yield),
                data.frame(curves[, c("sel", "F")], what = "SSB", value = curves$ssb))
  pts <- rbind(data.frame(msy[, c("sel", "F")], what = "Yield", value = msy$yield),
               data.frame(msy[, c("sel", "F")], what = "SSB", value = msy$ssb))
  ggplot2::ggplot(long, ggplot2::aes(.data$F, .data$value, colour = .data$sel)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(data = pts, size = 2.5) +
    ggplot2::facet_wrap(~what, scales = "free_y") +
    ggplot2::scale_colour_manual(values = floss_cols(length(sels))) +
    ggplot2::labs(x = "F (Fbar)", y = NULL, colour = "Selectivity",
                  title = sp$species,
                  subtitle = "Same stock, different selectivity; points: MSY") +
    ggplot2::theme_bw()
}

#' Spawning biomass versus vulnerable biomass by fleet
#'
#' Who reproduces and who is caught are not the same fish. The left panel
#' compares each fleet's selectivity-at-length with the maturity ogive; the
#' right panel compares the trend in spawning biomass (SSB, reproductive
#' potential) with the trend in each fleet's vulnerable biomass (the biomass
#' its catches and CPUE come from), both relative to the first year. With
#' dome-shaped selectivity, large spawners escape: a fleet's CPUE can fall
#' faster or slower than SSB, so an index can be a distorted picture of the
#' reproductive stock.
#'
#' @param om a `floss_om`.
#' @param iters iterations (medians are shown).
#' @return a `ggplot` (two panels, via facets).
#' @export
plot_ssb_vb <- function(om, iters = NULL) {
  sp <- om$species
  stk <- if (is.null(iters)) om$stock else FLCore::iter(om$stock, iters)
  yrs <- an(dimnames(stk)$year)
  med <- function(x) apply(matrix(c(x), nrow = length(yrs)), 1, stats::median)
  rel <- function(v) v / v[1]
  # trends: SSB and vulnerable biomass of each fleet (mid-year)
  tr <- data.frame(year = yrs, series = "SSB (spawners)", value = rel(med(ssb(stk))))
  for (f in om$fishery$fleets) {
    vb <- vuln_biomass(stk, c(om$sel[[f]]$sel_a), 0.5)
    for (p in fleet_periods(om, f)[-1]) {          # selectivity changes
      yy <- intersect(ac(yrs), ac(p$years))
      vb[, yy] <- vuln_biomass(stk, sel_a_for(om, p$sel), 0.5)[, yy]
    }
    tr <- rbind(tr, data.frame(year = yrs, series = paste("vulnerable to", f),
                               value = rel(med(vb))))
  }
  # selectivity at length vs maturity
  L <- seq(0, 1.2 * sp$lh$linf, length.out = 200)
  sel_len <- function(s) {
    v <- switch(s$type,
      logistic = FLRef::sel_logistic_len(L, L50 = s$L50, L95 = s$L95),
      normal = FLRef::sel_normal_len(L, Lpeak = s$Lpeak, sd = s$sd),
      dnormal = FLRef::sel_dnormal_len(L, Lpeak = s$Lpeak, sd_left = s$sd_left,
                                       sd_right = s$sd_right))
    v / max(v)
  }
  L95m <- sp$lh$lm * 1.15
  sl <- data.frame(x = L, series = "SSB (spawners)",
                   value = 1 / (1 + exp(-log(19) * (L - sp$lh$lm) / (L95m - sp$lh$lm))))
  for (f in om$fishery$fleets)
    sl <- rbind(sl, data.frame(x = L, series = paste("vulnerable to", f),
                               value = sel_len(sp$sel[[f]])))
  sl$panel <- "Maturity and selectivity at length (cm)"
  tr$panel <- "Trend relative to first year"
  names(tr)[1] <- "x"
  d <- rbind(sl, tr)
  d$panel <- factor(d$panel, levels = c("Maturity and selectivity at length (cm)",
                                        "Trend relative to first year"))
  cols <- c(`SSB (spawners)` = "#111111",
            stats::setNames(fleet_colours(om$fishery$fleets)[om$fishery$fleets],
                            paste("vulnerable to", om$fishery$fleets)))
  ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$value, colour = .data$series,
                                  linetype = .data$series == "SSB (spawners)")) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_linetype_manual(values = c(`TRUE` = 2, `FALSE` = 1), guide = "none") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL, colour = NULL, title = sp$species,
                  subtitle = "Who reproduces and who is caught: spawning vs vulnerable biomass") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
