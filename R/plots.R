# Plots for fisheries and operating models.

#' Plot fleet effort
#'
#' Median and 90% interval of nominal effort by fleet, with a few individual
#' iterations. With `efficiency = TRUE` the effective effort (nominal x
#' efficiency) is shown as well, to visualise effort creep.
#'
#' @param x a [fishery()] (effort is simulated) or a `floss_om`.
#' @param nsim iterations to simulate when `x` is a fishery.
#' @param iters iterations drawn as lines.
#' @param efficiency also show effective effort.
#' @return a `ggplot`.
#' @export
#' @examples
#' plot_effort(reef_fishery(creep = 0.02), efficiency = TRUE)
plot_effort <- function(x, nsim = 100, iters = 1:3, efficiency = FALSE) {
  eff <- if (inherits(x, "floss_om")) x[c("effort", "efficiency")] else
    effort_sim(x, nsim)
  df <- flqs_df(eff$effort, "Nominal")
  if (efficiency) {
    eeff <- FLQuants(Map(`*`, eff$effort, eff$efficiency))
    df <- rbind(df, flqs_df(eeff, "Effective"))
  }
  ribbon_plot(df, iters, "Relative effort") +
    ggplot2::facet_wrap(~qname, scales = "free_y", ncol = 1)
}

#' Plot apical fishing mortality or catch by fleet
#'
#' @param om a `floss_om`.
#' @param what `"fapic"` (apical F by fleet, relative to Fmsy) or `"catch"`
#'   (true catch by fleet, t).
#' @inheritParams plot_effort
#' @return a `ggplot`.
#' @export
plot_fleets <- function(om, what = c("fapic", "catch"), iters = 1:3) {
  what <- match.arg(what)
  if (what == "fapic") {
    fmsy <- c(om$stock@refpts["Fmsy"])
    df <- flqs_df(FLQuants(lapply(om$fapic, function(x) x / fmsy)), "")
    ylab <- expression(F[apical] / F[MSY])
  } else {
    df <- flqs_df(catch_fleet(om), "")
    ylab <- "Catch (t)"
  }
  ribbon_plot(df, iters, ylab) +
    ggplot2::facet_wrap(~qname, scales = "free_y") +
    ggplot2::ggtitle(om$species$species) +
    ggplot2::theme(legend.position = "none")
}

#' Plot stock status of an operating model
#'
#' Thin wrapper around [FLRef::plotAdvice()]. Several species are shown
#' relative to their own reference points (B/Bmsy, F/Fmsy, C/MSY).
#' @param om a `floss_om` or `floss_oms`.
#' @return a `ggplot`.
#' @export
plot_om <- function(om) {
  if (inherits(om, "floss_om")) return(FLRef::plotAdvice(om$stock) +
                                         ggplot2::ggtitle(om$species$species))
  stks <- FLStocks(lapply(om, `[[`, "stock"))
  names(stks) <- vapply(om, function(o) o$species$species, "")
  FLRef::plotAdvice(FLRef::stock2ratios(stks))
}

# helpers ---------------------------------------------------------------------

flqs_df <- function(x, type) {
  df <- as.data.frame(x)
  df$type <- type
  df
}

ribbon_plot <- function(df, iters, ylab) {
  q <- stats::aggregate(data ~ year + qname + type, df, function(v)
    stats::quantile(v, c(0.05, 0.5, 0.95)))
  q <- do.call(data.frame, q)
  names(q)[4:6] <- c("lo", "med", "hi")
  sub <- df[as.character(df$iter) %in% as.character(iters), ]
  multi <- length(unique(df$type)) > 1
  p <- ggplot2::ggplot(q, ggplot2::aes(.data$year)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi,
                                      fill = .data$type), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = .data$med, colour = .data$type),
                       linewidth = 0.9) +
    ggplot2::geom_line(data = sub, ggplot2::aes(
      y = .data$data, group = interaction(.data$iter, .data$type),
      colour = .data$type), linewidth = 0.3, alpha = 0.7) +
    ggplot2::labs(x = "Year", y = ylab, colour = NULL, fill = NULL) +
    ggplot2::theme_bw()
  if (!multi) p <- p + ggplot2::scale_colour_manual(values = "#0b5d8a") +
    ggplot2::scale_fill_manual(values = "#0b5d8a") +
    ggplot2::guides(colour = "none", fill = "none")
  else {
    cols <- stats::setNames(floss_cols(length(unique(df$type))), unique(df$type))
    p <- p + ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols)
  }
  p
}

#' @importFrom ggplot2 .data
NULL
