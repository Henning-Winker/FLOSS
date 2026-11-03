# Colours: the r4ss palette used in JABBA (rc4() and ss3col()), so that FLOSS
# plots look like JABBA and Stock Synthesis plots.

#' FLOSS colour palette (r4ss colours, as in JABBA)
#'
#' The "rich colours" palette of r4ss and JABBA (`rc4()`, `ss3col()`, after
#' `rich.colors` by Arni Magnusson): dark blue through green and yellow to
#' red. For three colours it returns blue, red and green3, as in JABBA.
#' Available for users who prefer the r4ss look; FLOSS plots use their own
#' muted palette by default.
#'
#' @param n number of colours.
#' @param alpha transparency (0-1).
#' @return a character vector of colours.
#' @export
#' @examples
#' floss.r4col(4)
#' barplot(rep(1, 6), col = floss.r4col(6))
floss.r4col <- function(n, alpha = 1) {
  n <- as.integer(n)
  if (n < 1) return(character(0))
  rich <- function(n) {
    x <- seq(0, 1, length.out = n)
    r <- 1 / (1 + exp(20 - 35 * x))
    g <- pmin(pmax(0, -0.8 + 6 * x - 5 * x^2), 1)
    b <- stats::dnorm(x, 0.25, 0.15) / max(stats::dnorm(x, 0.25, 0.15))
    grDevices::rgb(r, g, b)
  }
  col <- if (n > 3) rich(n + 1)[-1] else if (n < 3) rich(n) else c("blue", "red", "green3")
  if (alpha < 1) col <- grDevices::adjustcolor(col, alpha.f = alpha)
  col
}

# The FLOSS palette: muted, distinct colours for fleets, iterations and
# treatment runs; the survey is dark navy.
floss_pal <- c("#0A6E8C", "#A85400", "#5B45A8", "#2E7D32", "#8E2C48", "#B08900",
               "#4E6E81", "#C0587E")

# n colours from the FLOSS palette (beyond 8, interpolated)
floss_cols <- function(n) {
  if (n <= length(floss_pal)) floss_pal[seq_len(n)] else grDevices::colorRampPalette(floss_pal)(n)
}

# Fleet colours (fleets in fishery order, then the survey), shared by all plots
fleet_colours <- function(fleets) {
  c(stats::setNames(floss_cols(length(fleets)), fleets), survey = "#102330",
    `all fleets` = "#3A3A3A")
}

# Iteration colours, used when gears are already separated by panels
iter_colours <- function(iters) stats::setNames(floss_cols(length(iters)), iters)

# Run colours for treatments: baseline black, then the palette
run_colours <- function(runs) {
  stats::setNames(c("#111111", floss_cols(max(1, length(runs) - 1)))[seq_along(runs)], runs)
}

# Season colours: a sequential scale through the year
season_colours <- function(n) {
  grDevices::hcl.colors(n + 1, "viridis")[seq_len(n)]   # drop the palest yellow
}
