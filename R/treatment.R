# Treatments: controlled changes to an operating model, run against a
# baseline with common random numbers.

#' Define a treatment
#'
#' A treatment is one or more controlled changes applied to an operating model
#' from a given year. Run with [run_treatment()], the baseline and the treated
#' operating models share every random number (effort noise, recruitment
#' deviations, sampling), so all differences come from the treatment.
#'
#' Components (each optional, `NULL` = no change):
#'
#' * `effort`: list of `fleet`, `from`, `mult`: nominal effort of the fleet is
#'   multiplied by `mult` from year `from` (`mult = 0`: fleet leaves; a value
#'   > 0 for a fleet with no effort before has no effect, so use the fishery
#'   trend to add a fleet).
#' * `recruitment`: list with `type` = `"pulse"` (`years`, `mult`), `"trend"`
#'   (`from`, `rate` per year, e.g. -0.05) or `"regime"` (`from`, `mult`):
#'   multiplies the recruitment deviations.
#' * `selectivity`: list of `fleet`, `from`, `sel` (a [selectivity()]): the
#'   fleet's selectivity changes from year `from` (e.g. a gear switch or a
#'   minimum size).
#' * `creep`: list of `fleet`, `from`, `rate`: efficiency increases by `rate`
#'   per year from `from` (in addition to the fishery's own creep).
#' * `M`: list of `from`, `mult`: natural mortality is multiplied from `from`.
#'
#' Several entries of one component are given as a list of lists.
#'
#' @param name label.
#' @param effort,recruitment,selectivity,creep,M treatment components.
#' @return a list of class `floss_treatment`.
#' @export
#' @examples
#' treatment("halve linefish effort",
#'           effort = list(fleet = "linefish", from = 2010, mult = 0.5))
#' treatment("strong 2015 year class",
#'           recruitment = list(type = "pulse", years = 2015, mult = 5))
treatment <- function(name = "treatment", effort = NULL, recruitment = NULL,
                      selectivity = NULL, creep = NULL, M = NULL) {
  wrap <- function(x) if (is.null(x) || is.list(x[[1]])) x else list(x)
  tr <- structure(list(name = name, effort = wrap(effort),
                       recruitment = wrap(recruitment),
                       selectivity = wrap(selectivity), creep = wrap(creep),
                       M = M), class = "floss_treatment")
  for (s in tr$selectivity)
    if (!inherits(s$sel, "floss_sel")) stop("selectivity$sel must be a selectivity()")
  for (r in tr$recruitment)
    if (!r$type %in% c("pulse", "trend", "regime"))
      stop("recruitment type must be 'pulse', 'trend' or 'regime'")
  tr
}

#' @exportS3Method base::print
print.floss_treatment <- function(x, ...) {
  cat("FLOSS treatment:", x$name, "\n")
  for (e in x$effort) cat("  effort:", e$fleet, "x", e$mult, "from", e$from, "\n")
  for (r in x$recruitment) cat("  recruitment:", r$type,
    switch(r$type, pulse = paste("x", r$mult, "in", paste(r$years, collapse = ",")),
           trend = paste(r$rate * 100, "% per year from", r$from),
           regime = paste("x", r$mult, "from", r$from)), "\n")
  for (s in x$selectivity) cat("  selectivity:", s$fleet, "->", s$sel$type,
                               "from", s$from, "\n")
  for (c in x$creep) cat("  creep:", c$fleet, c$rate * 100, "% per year from",
                         c$from, "\n")
  if (!is.null(x$M)) cat("  M x", x$M$mult, "from", x$M$from, "\n")
  invisible(x)
}

#' Run a baseline and treatments with common random numbers
#'
#' Builds the baseline operating model(s) and one per treatment, all with the
#' same effort noise, recruitment deviations and (when `oem` is given)
#' sampling random numbers.
#'
#' @param species species ids or [species()] objects.
#' @param fishery a [fishery()].
#' @param treatments a [treatment()] or a named list of treatments.
#' @param nsim iterations.
#' @param oem optional [oem_settings()]: also sample data for `iters`.
#' @param iters iterations to sample.
#' @param seed base seed (operating models and sampling).
#' @param ... passed to [build_om()].
#' @return a list of class `floss_treatments` with `oms` and (if `oem`)
#'   `data`, each a named list (`baseline`, then the treatments).
#' @export
run_treatment <- function(species, fishery, treatments, nsim = 20, oem = NULL,
                          iters = 1, seed = fishery$seed, ...) {
  if (inherits(treatments, "floss_treatment"))
    treatments <- stats::setNames(list(treatments), treatments$name)
  if (is.null(names(treatments)))
    names(treatments) <- vapply(treatments, `[[`, "", "name")
  sps <- lapply(if (inherits(species, "floss_species")) list(species) else species,
                function(s) if (is.character(s)) load_species(s) else s)
  effort <- effort_sim(fishery, nsim)
  one <- function(tr) {
    oms <- lapply(sps, function(s)
      build_om(s, fishery, nsim = nsim, effort = effort, seed = seed,
               treatment = tr, ...))
    names(oms) <- vapply(sps, `[[`, "", "id")
    structure(oms, class = "floss_oms")
  }
  oms <- c(list(baseline = one(NULL)), lapply(treatments, one))
  out <- list(oms = oms, treatments = treatments)
  if (!is.null(oem))
    out$data <- lapply(oms, sample_data, oem = oem, iters = iters, seed = seed)
  structure(out, class = "floss_treatments")
}

#' @exportS3Method base::print
print.floss_treatments <- function(x, ...) {
  cat("FLOSS treatment comparison:", length(x$oms) - 1, "treatment(s) vs baseline\n")
  for (n in names(x$treatments)) print(x$treatments[[n]])
  if (!is.null(x$data)) cat("data sampled for iterations",
                            paste(x$data[[1]]$settings$iters, collapse = ", "), "\n")
  invisible(x)
}

#' Plot the response to treatments
#'
#' True B/Bmsy, F/Fmsy, catch/MSY and recruitment (relative to the baseline
#' mean) for the baseline and each treatment: median and 90% interval over
#' iterations; the treatment start is marked.
#' @param x a `floss_treatments`.
#' @param species species id; default the first.
#' @return a `ggplot`.
#' @export
plot_treatment <- function(x, species = NULL) {
  species <- species %||% names(x$oms$baseline)[1]
  base_rec <- mean(c(rec(x$oms$baseline[[species]]$stock)))
  d <- do.call(rbind, lapply(names(x$oms), function(n) {
    o <- x$oms[[n]][[species]]
    s <- o$stock
    rp <- x$oms$baseline[[species]]$stock@refpts
    q <- function(v, nm) {
      m <- matrix(c(v), nrow = dims(s)$year)
      data.frame(year = an(dimnames(s)$year), qname = nm,
                 lo = apply(m, 1, stats::quantile, 0.05),
                 med = apply(m, 1, stats::median),
                 hi = apply(m, 1, stats::quantile, 0.95))
    }
    rbind(q(ssb(s) / c(rp["Bmsy"]), "B/Bmsy"), q(fbar(s) / c(rp["Fmsy"]), "F/Fmsy"),
          q(catch(s) / c(rp["MSY"]), "Catch/MSY"), q(rec(s) / base_rec, "Recruitment")) |>
      transform(run = n)
  }))
  d$qname <- factor(d$qname, levels = c("B/Bmsy", "F/Fmsy", "Catch/MSY", "Recruitment"))
  d$run <- factor(d$run, levels = names(x$oms))
  starts <- unique(unlist(lapply(x$treatments, function(t)
    c(vapply(t$effort, `[[`, 0, "from"), vapply(t$selectivity, `[[`, 0, "from"),
      vapply(t$creep, `[[`, 0, "from"), t$M$from,
      unlist(lapply(t$recruitment, function(r) r$from %||% min(r$years)))))))
  cols <- run_colours(levels(d$run))   # baseline always black
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$med, colour = .data$run,
                                  fill = .data$run)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = "grey40") +
    ggplot2::geom_vline(xintercept = starts, linetype = 2, colour = "grey50") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi),
                         colour = NA, alpha = 0.15) +
    ggplot2::geom_line(data = d[d$run != "baseline", ], linewidth = 0.9) +
    # baseline on top, dashed, so treatments that leave a quantity unchanged
    # do not hide it
    ggplot2::geom_line(data = d[d$run == "baseline", ], linewidth = 0.8,
                       linetype = 2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::facet_wrap(~qname, scales = "free_y") +
    ggplot2::labs(x = "Year", y = NULL, colour = NULL, fill = NULL,
                  title = x$oms$baseline[[species]]$species$species,
                  subtitle = "Median and 90% interval; dashed: treatment starts; same random numbers in all runs") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

# ---- one-gear teaching scenarios ---------------------------------------------

#' One-gear fishery with constant (or stepped) F/Fmsy
#'
#' A fishery with a single fleet whose effort is directly F/Fmsy when used
#' with [one_gear_species()]. Use for light vs hard fishing comparisons.
#'
#' @param F_rel F/Fmsy: one value (constant after a short ramp), or a vector
#'   with one value per year.
#' @param years years.
#' @param fleet fleet name.
#' @param ramp years to ramp up from near zero to `F_rel` (0 = immediately).
#' @param sigmaE effort noise.
#' @param seed seed.
#' @return a `floss_fishery`.
#' @export
#' @examples
#' light <- one_gear_fishery(0.5)
#' hard  <- one_gear_fishery(2)
one_gear_fishery <- function(F_rel = 1, years = 1969:2025, fleet = "handline",
                             ramp = 10, sigmaE = 0.1, seed = 123) {
  ny <- length(years)
  v <- if (length(F_rel) == ny) F_rel else {
    r <- min(ramp, ny - 1)
    c(if (r > 0) seq(0.05, F_rel, length.out = r + 1)[-(r + 1)], rep(F_rel, ny - r))
  }
  fishery(fleets = fleet, years = years,
          trends = stats::setNames(list(effort_trend("values", values = v)), fleet),
          sigmaE = sigmaE, name = paste0("one_gear_", fleet), seed = seed)
}

#' Adapt a species to a one-gear fishery
#'
#' Copies a species with a single fleet (default a logistic handline with
#' L50 at the length at maturity) and scales its catchability so that the
#' fishery's effort equals F/Fmsy.
#' @param sp a [species()].
#' @param fleet fleet name (as in [one_gear_fishery()]).
#' @param sel a [selectivity()]; default logistic with L50 = `lm`, L95 = 1.2
#'   `lm`.
#' @return a `floss_species`.
#' @export
one_gear_species <- function(sp, fleet = "handline", sel = NULL) {
  if (is.null(sel)) sel <- selectivity("logistic", L50 = sp$lh$lm, L95 = 1.2 * sp$lh$lm)
  s1 <- species(sp$species, sp$common_name, lh = sp$lh, msy = sp$msy,
                sel = stats::setNames(list(sel), fleet),
                relq = stats::setNames(1, fleet), productivity = sp$productivity,
                id = sp$id)
  # apical F = relq x Fmsy(fbar) x effort; fbar = apical x mean selectivity over
  # the fbar ages, so relq = 1 / mean(sel) makes effort = F/Fmsy
  sa <- c(species_sel(s1)[[fleet]]$sel_a)
  fb <- seq(1, sp$lh$tmax - 1) + 1          # ages 1..tmax-1 (index from age 0)
  s1$relq[[fleet]] <- 1 / mean(sa[fb])
  s1
}

# ---- applying a treatment inside build_om ------------------------------------

apply_effort_treatment <- function(effort, tr, years) {
  if (is.null(tr)) return(effort)
  for (e in tr$effort) {
    yy <- ac(years[years >= e$from])
    effort$effort[[e$fleet]][, yy] <- effort$effort[[e$fleet]][, yy] * e$mult
  }
  for (c in tr$creep) {
    yy <- years[years >= c$from]
    eta <- effort$efficiency[[c$fleet]]
    eta[, ac(yy)] <- eta[, ac(yy)] * rep((1 + c$rate)^(yy - c$from), dims(eta)$iter)
    effort$efficiency[[c$fleet]] <- eta
  }
  effort
}

rec_multiplier <- function(tr, years) {
  m <- rep(1, length(years))
  if (is.null(tr)) return(m)
  for (r in tr$recruitment) {
    if (r$type == "pulse") m[years %in% r$years] <- m[years %in% r$years] * r$mult
    if (r$type == "regime") m[years >= r$from] <- m[years >= r$from] * r$mult
    if (r$type == "trend") {
      k <- years >= r$from
      m[k] <- m[k] * exp(r$rate * (years[k] - r$from))
    }
  }
  m
}

# Selectivity spec of a fleet in each year: list of periods (years, sel)
fleet_periods <- function(om, fleet) {
  yrs <- an(dimnames(om$stock)$year)
  chg <- Filter(function(s) identical(s$fleet, fleet), om$treatment$selectivity)
  if (!length(chg)) return(list(list(years = yrs, sel = om$species$sel[[fleet]])))
  from <- chg[[1]]$from
  out <- list(list(years = yrs[yrs < from], sel = om$species$sel[[fleet]]),
              list(years = yrs[yrs >= from], sel = chg[[1]]$sel))
  Filter(function(p) length(p$years) > 0, out)
}
