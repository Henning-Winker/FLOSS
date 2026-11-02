# Fleet effort: deterministic trends, shared process noise and effort creep.

#' Deterministic effort trend values
#'
#' @param trend a [effort_trend()].
#' @param years integer vector of years.
#' @return numeric vector, one value per year.
#' @export
#' @examples
#' trend_values(effort_trend(fhi = 2, flo = 1, breaks = c(0.5, 0.8)), 1969:2025)
trend_values <- function(trend, years) {
  x <- as.numeric(years)
  if (trend$type == "values") {
    if (length(trend$values) != length(x)) stop("trend has wrong length")
    return(trend$values)
  }
  # up-down-constant shape, as FLRef::fudc() but with an explicit start e0
  steps <- length(x)
  y0 <- x[1]
  y1 <- x[max(2, floor(steps * trend$breaks[1]))]
  y2 <- x[min(steps, max(ceiling(steps * trend$breaks[2]),
                         match(y1, x) + 1))]
  hi <- trend$fhi * trend$fref
  lo <- trend$flo * trend$fref
  out <- rep(trend$e0, steps)
  up <- x > y0 & x <= y1
  dn <- x > y1 & x <= y2
  out[up] <- (hi - trend$e0) / (y1 - y0) * (x[up] - y0) + trend$e0
  out[dn] <- hi - (hi - lo) / (y2 - y1) * (x[dn] - y1)
  out[x > y2] <- lo
  out
}

#' Simulate fleet effort
#'
#' Nominal effort \eqn{E_{f,y}} is the trend times multiplicative lognormal
#' noise with mean one, \eqn{\exp(\epsilon - \sigma_E^2/2)},
#' \eqn{\epsilon\sim N(0,\sigma_E^2)}. The noise is drawn once per fishery and
#' iteration, so all species share it. Efficiency (effort creep) is
#' \eqn{\eta_{f,y} = (1 + c_f)^{\max(0, y - y^c_f)}}; it scales fishing
#' mortality but not the nominal effort index.
#'
#' @param fishery a [fishery()].
#' @param nsim number of iterations.
#' @param seed random seed; default the fishery's seed, so repeated calls
#'   give the same effort.
#' @return a list with `FLQuants` `effort` (nominal) and `efficiency`, one
#'   element per fleet, each `FLQuant(age = "all", year, iter = nsim)`.
#' @export
#' @examples
#' eff <- effort_sim(reef_fishery(), nsim = 5)
#' eff$effort$traps
effort_sim <- function(fishery, nsim = 100, seed = fishery$seed) {
  fishery <- validate_fishery(fishery)
  yrs <- fishery_years(fishery)
  ny <- length(yrs)
  with_seed(seed, {
    eff <- lapply(fishery$fleets, function(f) {
      sig <- fishery$sigmaE[[f]]
      dev <- matrix(stats::rlnorm(ny * nsim, -0.5 * sig^2, sig), ny, nsim)
      E <- flq(trend_values(fishery$trends[[f]], yrs) * dev, years = yrs,
               nsim = nsim, units = "relative")
      eta <- (1 + fishery$creep[[f]])^pmax(0, yrs - fishery$creep_start[[f]])
      list(E = E, eta = flq(rep(eta, nsim), years = yrs, nsim = nsim,
                            units = "relative"))
    })
  })
  names(eff) <- fishery$fleets
  list(effort = FLQuants(lapply(eff, `[[`, "E")),
       efficiency = FLQuants(lapply(eff, `[[`, "eta")))
}

# Evaluate expr with a fixed seed and restore the caller's RNG state.
with_seed <- function(seed, expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

# A reproducible integer seed for a (base seed, tag) pair, e.g. one per species.
seed_for <- function(seed, tag) {
  h <- sum(utf8ToInt(tag) * seq_along(utf8ToInt(tag))) %% 100000L
  as.integer((as.numeric(seed) + h) %% .Machine$integer.max)
}

# Build an FLQuant(age, year, iter) from a vector ordered age, year, iter.
flq <- function(x, age = "all", years, nsim, units = "NA") {
  dn <- list(age = as.character(age), year = as.character(years),
             unit = "unique", season = "all", area = "unique",
             iter = as.character(seq_len(nsim)))
  FLQuant(array(as.numeric(x), dim = lengths(dn), dimnames = dn), units = units)
}
