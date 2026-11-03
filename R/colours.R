# Colours: the r4ss palette used in JABBA (rc4() and ss3col()), so that FLOSS
# plots look like JABBA and Stock Synthesis plots.

#' FLOSS colour palette (r4ss colours, as in JABBA)
#'
#' The "rich colours" palette of r4ss and JABBA (`rc4()`, `ss3col()`, after
#' `rich.colors` by Arni Magnusson): dark blue through green and yellow to
#' red. For three colours it returns blue, red and green3, as in JABBA.
#' FLOSS uses it for fleets and survey, seasons, selectivities and treatment
#' runs.
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

# Fleet colours (fleets first, then the survey), shared by all data plots
fleet_colours <- function(fleets) {
  cols <- floss.r4col(length(fleets) + 1)
  c(stats::setNames(cols[seq_along(fleets)], fleets), survey = cols[length(fleets) + 1])
}

# Run colours for treatments: baseline black, then the palette
run_colours <- function(runs) {
  stats::setNames(c("#111111", floss.r4col(max(1, length(runs) - 1)))[seq_along(runs)], runs)
}
