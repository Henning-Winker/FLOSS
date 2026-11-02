# Fishery specification: fleets, years and fleet-level nominal effort.
# The fishery is species-independent; all species share one effort realisation.

#' Effort trend for one fleet
#'
#' Describes the shape of a fleet's nominal effort over time. `"fudc"` is the
#' up-down-constant shape of [FLRef::fudc()]: effort starts at `e0` in the
#' first year, rises linearly to `fhi x fref` by the first break, moves
#' linearly to `flo x fref` by the second break and stays there. `"values"` takes a user-supplied vector, one value per year.
#'
#' @param type `"fudc"` or `"values"`.
#' @param e0 effort in the first year (relative units; near-unfished start).
#' @param fhi,flo,fref levels for `"fudc"` (relative units).
#' @param breaks two proportions of the time series (0-1) where the trend
#'   reaches `fhi` and then `flo`.
#' @param values numeric vector of relative effort, one per year (`"values"`).
#' @return a list of class `floss_trend`.
#' @export
#' @examples
#' effort_trend(fhi = 1, flo = 0.8, breaks = c(0.4, 0.6))
effort_trend <- function(type = c("fudc", "values"), e0 = 0.05, fhi = 1,
                         flo = 0.8, fref = 1, breaks = c(0.5, 0.75),
                         values = NULL) {
  type <- match.arg(type)
  out <- switch(type,
    fudc = list(type = "fudc", e0 = e0, fhi = fhi, flo = flo, fref = fref,
                breaks = as.numeric(breaks)),
    values = {
      if (is.null(values)) stop("type = 'values' needs 'values'")
      list(type = "values", values = as.numeric(values))
    })
  structure(out, class = "floss_trend")
}

#' Define a fishery
#'
#' A fishery is a set of fleets that share a time horizon. Each fleet has a
#' nominal effort trend, process noise on effort that is **common to all
#' species** (so species caught by the same fleet share effort anomalies) and
#' an optional annual efficiency increase (effort creep) that raises fishing
#' mortality but is not visible in the nominal effort index.
#'
#' Fishing mortality of species \eqn{s} by fleet \eqn{f} in year \eqn{y} is
#' \deqn{F_{s,f,y} = q_{s,f}\,\eta_{f,y}\,E_{f,y},\quad
#'       \eta_{f,y} = (1 + c_f)^{\max(0,\,y - y^{c}_f)},}
#' where \eqn{E} is nominal effort (trend x lognormal noise), \eqn{c_f} the
#' creep rate and \eqn{q_{s,f}} the species catchability, set in
#' [species()] as relative catchability x \eqn{F_{MSY}}.
#'
#' @param fleets character vector of fleet names.
#' @param years integer vector of years (first year = unfished start).
#' @param trends named list of [effort_trend()] objects, one per fleet.
#' @param sigmaE lognormal SD of the effort process noise (shared across
#'   species); scalar or one value per fleet.
#' @param creep annual efficiency increase per fleet (0.02 = 2% per year).
#' @param creep_start first year of creep per fleet (default: first year).
#' @param name a label.
#' @param seed seed for the effort realisation; every species built with this
#'   fishery sees the same effort.
#' @return a list of class `floss_fishery`.
#' @export
#' @examples
#' fishery(fleets = c("traps", "linefish"), years = 1980:2025,
#'         trends = list(traps = effort_trend(fhi = 1, flo = 0.8),
#'                       linefish = effort_trend(fhi = 2, flo = 1)))
fishery <- function(fleets, years, trends, sigmaE = 0.15, creep = 0,
                    creep_start = NULL, name = "fishery", seed = 123) {
  fleets <- as.character(fleets)
  nf <- length(fleets)
  years <- as.integer(years)
  if (is.null(names(trends))) names(trends) <- fleets
  per_fleet <- function(x, what) {
    if (length(x) == 1) x <- rep(x, nf)
    if (length(x) != nf) stop("'", what, "' needs 1 or ", nf, " values")
    if (is.null(names(x))) names(x) <- fleets
    as.list(x[fleets])
  }
  if (is.null(creep_start)) creep_start <- min(years)
  trends <- lapply(trends[fleets], function(t) {
    if (!inherits(t, "floss_trend")) t <- do.call(effort_trend, t)
    t
  })
  out <- structure(list(
    name = name, version = floss_spec_version(),
    fleets = fleets, years = range(years),
    trends = trends,
    sigmaE = per_fleet(sigmaE, "sigmaE"),
    creep = per_fleet(creep, "creep"),
    creep_start = per_fleet(as.integer(creep_start), "creep_start"),
    seed = as.integer(seed)), class = "floss_fishery")
  validate_fishery(out)
}

#' @rdname fishery
#' @param x a `floss_fishery`.
#' @export
validate_fishery <- function(x) {
  stopifnot(inherits(x, "floss_fishery"))
  if (anyDuplicated(x$fleets)) stop("fleet names must be unique")
  miss <- setdiff(x$fleets, names(x$trends))
  if (length(miss)) stop("no effort trend for fleet(s): ", toString(miss))
  if (diff(x$years) < 3) stop("need at least 4 years")
  nyr <- diff(x$years) + 1
  for (f in x$fleets) {
    t <- x$trends[[f]]
    if (t$type == "values" && length(t$values) != nyr)
      stop("trend for '", f, "' has ", length(t$values), " values, need ", nyr)
    if (t$type == "fudc" && (length(t$breaks) != 2 || any(t$breaks <= 0) ||
                             any(t$breaks >= 1) || t$breaks[1] > t$breaks[2]))
      stop("breaks for '", f, "' must be two increasing proportions in (0, 1)")
  }
  if (any(unlist(x$sigmaE) < 0)) stop("sigmaE must be >= 0")
  x
}

#' @rdname fishery
#' @export
fishery_years <- function(x) seq(x$years[1], x$years[2])

#' The three-fleet reef fishery used in the examples
#'
#' Traps, gillnet and linefish, 1969-2025, with the effort trends of the
#' multi-gear reef operating model (Winker, FLsim multigear reef example).
#'
#' @inheritParams fishery
#' @return a `floss_fishery`.
#' @export
#' @examples
#' reef_fishery()
reef_fishery <- function(sigmaE = 0.15, creep = 0, seed = 123) {
  fishery(
    name = "reef_3fleets",
    fleets = c("traps", "gillnet", "linefish"),
    years = 1969:2025,
    trends = list(
      traps    = effort_trend(fhi = 1.0, flo = 0.8, breaks = c(0.4, 0.6)),
      gillnet  = effort_trend(fhi = 0.5, flo = 1.5, breaks = c(0.5, 0.8)),
      linefish = effort_trend(fhi = 2.5, flo = 0.9, breaks = c(0.6, 0.8))),
    sigmaE = sigmaE, creep = creep, seed = seed)
}

#' @exportS3Method base::print
print.floss_fishery <- function(x, ...) {
  cat("FLOSS fishery '", x$name, "': ", length(x$fleets), " fleets, ",
      x$years[1], "-", x$years[2], "\n", sep = "")
  tab <- data.frame(
    fleet = x$fleets,
    trend = vapply(x$trends, function(t) if (t$type == "fudc")
      sprintf("fudc(e0=%g, hi=%g, lo=%g, breaks=%s)", t$e0, t$fhi, t$flo,
              paste(t$breaks, collapse = "/")) else "values", ""),
    sigmaE = unlist(x$sigmaE), creep = unlist(x$creep),
    creep_start = unlist(x$creep_start), row.names = NULL)
  print(tab, row.names = FALSE)
  invisible(x)
}

floss_spec_version <- function() "0.1"
