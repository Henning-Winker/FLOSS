# Seasonal sampling: length and biological samples are taken in seasons
# within the year. Fish grow through the year, and the catch-at-age of a
# season reflects mortality up to that season, so the length samples of
# different seasons differ. When the sample is spread unevenly over seasons
# (imbalanced sampling), the pooled annual length distribution is biased.

#' Seasonal sampling design
#'
#' Describes how the annual length and biological samples of a fleet are
#' spread over seasons. The default, four seasons with equal shares, is a
#' balanced design. Imbalanced designs:
#'
#' * **fixed imbalance**: unequal `share`, e.g. `c(0.5, 0.3, 0.2, 0)` (most
#'   samples early in the year, none in season 4, e.g. a monsoon closure of
#'   the landing site);
#' * **shift**: the shares change from a given year (`shift`), e.g. a new
#'   observer programme that samples in the second half of the year only;
#' * **opportunistic**: the shares vary at random from year to year around
#'   `share` (`kappa`, the concentration of a Dirichlet distribution: 50 =
#'   small departures, 5 = large ones), as when samples are taken whenever
#'   an observer happens to be available.
#'
#' Within a season, the age composition of the catch is the Baranov catch
#' between the start and end of the season (fishing and natural mortality
#' constant through the year), and lengths are those at the middle of the
#' season (von Bertalanffy growth with age + season mid-point).
#'
#' @param n number of seasons (1 = one annual sample at mid-year, as in
#'   FLOSS <= 0.4).
#' @param share sample share by season (rescaled to sum to 1); default equal.
#' @param shift optional list of `from` (year) and `share`: the shares from
#'   that year on. Several shifts as a list of such lists.
#' @param kappa Dirichlet concentration for random year-to-year shares; `Inf`
#'   (default) = no random variation.
#' @param labels season names; default `"S1"`, `"S2"`, ...
#' @return a list of class `floss_seasons`.
#' @export
#' @examples
#' season_design()                                   # balanced, 4 seasons
#' season_design(share = c(0.6, 0.3, 0.1, 0))         # fixed imbalance
#' season_design(shift = list(from = 2015, share = c(0, 0, 0.5, 0.5)))
#' season_design(kappa = 5)                           # opportunistic
season_design <- function(n = 4, share = NULL, shift = NULL, kappa = Inf,
                          labels = NULL) {
  n <- as.integer(n)
  if (n < 1 || n > 12) stop("n must be between 1 and 12")
  norm <- function(s, what) {
    s <- as.numeric(s)
    if (length(s) != n) stop(what, " needs ", n, " values (one per season)")
    if (any(s < 0) || sum(s) <= 0) stop(what, " must be >= 0 with a positive sum")
    s / sum(s)
  }
  share <- if (is.null(share)) rep(1 / n, n) else norm(share, "share")
  if (!is.null(shift) && !is.null(shift$from)) shift <- list(shift)
  shift <- lapply(shift, function(s) list(from = as.integer(s$from), share = norm(s$share, "shift share")))
  if (length(shift)) shift <- shift[order(vapply(shift, `[[`, 0L, "from"))]
  if (!is.numeric(kappa) || kappa <= 0) stop("kappa must be positive (Inf = no random variation)")
  structure(list(n = n, share = share, shift = shift, kappa = kappa,
                 labels = labels %||% paste0("S", seq_len(n))),
            class = "floss_seasons")
}

#' @exportS3Method base::print
print.floss_seasons <- function(x, ...) {
  cat("FLOSS seasonal sampling:", x$n, "season(s)",
      if (is_balanced(x)) "(balanced)" else "(imbalanced)", "\n")
  cat("  share:", paste(sprintf("%.2f", x$share), collapse = " "), "\n")
  for (s in x$shift) cat("  from", s$from, ":", paste(sprintf("%.2f", s$share), collapse = " "), "\n")
  if (is.finite(x$kappa)) cat("  random year-to-year shares, kappa =", x$kappa, "\n")
  invisible(x)
}

is_balanced <- function(x) {
  eq <- function(s) isTRUE(all.equal(s, rep(1 / x$n, x$n)))
  eq(x$share) && all(vapply(x$shift, function(s) eq(s$share), TRUE)) && !is.finite(x$kappa)
}

# Mid-point of each season as a fraction of the year
season_mid <- function(n) (seq_len(n) - 0.5) / n

# Shares by year (rows) and season (columns); random part uses the current RNG
season_shares <- function(design, years) {
  years <- as.integer(years)
  m <- t(vapply(years, function(y) {
    s <- design$share
    for (p in design$shift) if (y >= p$from) s <- p$share
    s
  }, numeric(design$n)))
  if (design$n == 1) m <- matrix(1, length(years), 1)
  if (is.finite(design$kappa) && design$n > 1) {
    m <- t(apply(m, 1, function(s) {
      g <- ifelse(s > 0, stats::rgamma(length(s), shape = design$kappa * s), 0)
      if (sum(g) > 0) g / sum(g) else s
    }))
  }
  dimnames(m) <- list(year = years, season = design$labels)
  m
}

# Split a sample size over seasons (largest remainder, sums exactly to n)
split_n <- function(n, share) {
  raw <- n * share
  k <- floor(raw)
  rest <- n - sum(k)
  if (rest > 0) k[order(raw - k, decreasing = TRUE)[seq_len(rest)]] <-
    k[order(raw - k, decreasing = TRUE)[seq_len(rest)]] + 1
  k
}

# Share of the annual catch-at-age taken in each season (age x season), for a
# year and iteration: Baranov catch between season start and end
season_catch_frac <- function(z, n) {
  t0 <- (seq_len(n) - 1) / n
  t1 <- seq_len(n) / n
  tot <- 1 - exp(-z)
  vapply(seq_len(n), function(s) ifelse(tot > 0, (exp(-z * t0[s]) - exp(-z * t1[s])) / tot, 1 / n),
         numeric(length(z)))
}

# Conditional age-length key of a gear at a time of year
cond_alk_at <- function(gear, timing) {
  sl <- as.numeric(gear$sel_len)
  names(sl) <- dimnames(gear$sel_len)$len
  ialk <- FLRef::gear_ialk(gear, age = as.numeric(rownames(gear$condALK)), timing = timing)
  FLRef::condition_alk(ialk, sl)
}

# Seasonal design of a fleet (one design for all fleets, or a list by fleet)
fleet_seasons <- function(oem, fleet) {
  s <- oem$seasons
  if (is.null(s)) return(season_design(1))
  if (inherits(s, "floss_seasons")) return(s)
  s[[fleet]] %||% season_design(1)
}

#' Pool seasonal length samples into annual ones
#'
#' Sums the length samples of all seasons of a year. With `weight = "equal"`
#' each season gets the same weight whatever its sample size (each season's
#' distribution is rescaled to the mean seasonal sample size before
#' summing), which removes the bias of imbalanced sampling when the catch
#' is spread evenly over the year.
#'
#' @param lfd the `lfd` table of a `floss_data` (with or without `season`).
#' @param weight `"none"` (sum the counts) or `"equal"` (equal season weights).
#' @return the table without the `season` column.
#' @export
lfd_annual <- function(lfd, weight = c("none", "equal")) {
  weight <- match.arg(weight)
  if (is.null(lfd) || !nrow(lfd) || is.null(lfd$season)) return(lfd)
  key <- c("species", "iter", "year", "source")
  if (weight == "equal") {
    g <- interaction(lfd[, c(key, "season")], drop = TRUE)
    tot_s <- stats::ave(lfd$n, g, FUN = sum)
    gy <- interaction(lfd[, key], drop = TRUE)
    ns <- stats::ave(as.numeric(lfd$season), gy, FUN = function(v) length(unique(v)))
    tot_y <- stats::ave(lfd$n, gy, FUN = sum)
    lfd$n <- lfd$n / tot_s * tot_y / ns
  }
  out <- stats::aggregate(n ~ species + iter + year + source + len_mm, lfd, sum)
  if (weight == "none") out$n <- as.integer(round(out$n))
  out[order(out$species, out$iter, out$year, out$source, out$len_mm), ]
}
