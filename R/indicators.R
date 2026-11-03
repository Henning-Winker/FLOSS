# Size-based indicators and an indicator dashboard. These describe the data
# (and compare them with the truth); they are not assessment models.

#' Length-based indicators by year and source
#'
#' From the length samples of each year and source:
#'
#' * `Lc`: length of first full capture, here the mid-point of the modal
#'   length bin of the pooled samples of the source (or given).
#' * `Lmean`: mean length of fish at or above `Lc`; `L95`: 95th percentile;
#'   `Lmax5`: mean length of the largest 5% of fish.
#' * `Pmat`, `Popt`, `Pmega` (Froese 2004): proportion of fish above the
#'   length at maturity, within 10% of the optimum length
#'   \eqn{L_{opt} = 3L_\infty/(3 + M/k)}, and above 1.1 \eqn{L_{opt}}.
#' * `Z_lbar`: total mortality from mean length (Beverton and Holt 1956),
#'   \eqn{Z = k(L_\infty - \bar L)/(\bar L - L_c)}. It assumes equilibrium:
#'   constant recruitment and mortality, and flat-topped selectivity.
#' * `Z_true` (when the answer key is present): true Fbar + M in that year,
#'   to show when and why `Z_lbar` departs from it.
#'
#' @param x a `floss_data`, a [read_dataset()] list, or a `floss_treatments`
#'   (then one set of indicators per run, column `run`).
#' @param species species name.
#' @param iter iteration(s).
#' @param source sources; default all.
#' @param Lc length of first capture (cm) per source (named) or one value;
#'   default the modal length.
#' @param season_weight how seasonal length samples are pooled into a year
#'   ([lfd_annual()]): `"none"` sums the counts (imbalanced sampling then
#'   biases the indicators), `"equal"` gives every season the same weight.
#' @param lh optional list with `linf`, `k`, `M`, `lm` (cm, per year);
#'   default from the species settings or `lhpars`.
#' @return a data frame, one row per species, iteration, year and source.
#' @export
length_indicators <- function(x, species, iter = 1, source = NULL, Lc = NULL,
                              lh = NULL, season_weight = c("none", "equal")) {
  season_weight <- match.arg(season_weight)
  if (inherits(x, "floss_treatments")) {
    parts <- lapply(names(x$data), function(n)
      length_indicators(x$data[[n]], species, iter, source, Lc, lh, season_weight))
    out <- do.call(rbind, Map(function(p, n) cbind(run = n, p), parts, names(x$data)))
    out$run <- factor(out$run, levels = names(x$data))
    attr(out, "Lopt") <- attr(parts[[1]], "Lopt")   # cbind/rbind drop attributes
    attr(out, "lm") <- attr(parts[[1]], "lm")
    return(out)
  }
  lh <- lh %||% true_lh(x, species)
  if (is.null(lh)) stop("life-history parameters needed: pass 'lh'")
  l <- pick(lfd_annual(x$lfd, season_weight), species, iter, source)
  bw <- lfd_bin_mm(x, species) / 10
  l$L <- l$len_mm / 10 + bw / 2
  Lopt <- 3 * lh$linf / (3 + lh$M / lh$k)
  srcs <- unique(l$source)
  lc_src <- vapply(srcs, function(s) {
    if (!is.null(Lc)) return(if (length(Lc) == 1) Lc else Lc[[s]])
    ls <- stats::aggregate(n ~ L, l[l$source == s, ], sum)
    ls$L[which.max(ls$n)]
  }, 0)
  wq <- function(v, w, p) {                      # weighted quantile
    o <- order(v); v <- v[o]; cw <- cumsum(w[o]) / sum(w)
    v[which(cw >= p)[1]]
  }
  rows <- lapply(split(l, list(l$iter, l$year, l$source), drop = TRUE), function(d) {
    lc <- lc_src[[d$source[1]]]
    a <- d[d$L >= lc, ]
    lbar <- sum(a$L * a$n) / sum(a$n)
    l95 <- wq(d$L, d$n, 0.95)
    data.frame(species = species, iter = d$iter[1], year = d$year[1],
               source = d$source[1], n = sum(d$n), Lc = lc, Lmean = lbar,
               L95 = l95, Lmax5 = sum((d$L * d$n)[d$L >= l95]) / sum(d$n[d$L >= l95]),
               Pmat = sum(d$n[d$L >= lh$lm]) / sum(d$n),
               Popt = sum(d$n[d$L >= 0.9 * Lopt & d$L <= 1.1 * Lopt]) / sum(d$n),
               Pmega = sum(d$n[d$L > 1.1 * Lopt]) / sum(d$n),
               Z_lbar = lh$k * (lh$linf - lbar) / (lbar - lc))
  })
  out <- do.call(rbind, rows)
  if (!is.null(x$truth)) {
    tr <- x$truth[x$truth$species == species, ]
    out$Z_true <- tr$fbar[match(paste(out$iter, out$year), paste(tr$iter, tr$year))] + lh$M
  }
  out <- out[order(out$source, out$iter, out$year), ]
  rownames(out) <- NULL
  attr(out, "Lopt") <- Lopt
  attr(out, "lm") <- lh$lm
  out
}

#' Plot length-based indicators
#'
#' Mean length (with the optimum length and `Lc`), the share of mature,
#' optimum-sized and mega-spawning fish, and Z from mean length against the
#' true Z. With several runs (from a `floss_treatments`), runs are coloured
#' and sources are panels.
#'
#' @param li result of [length_indicators()].
#' @param source sources to show; default all.
#' @param stat with several iterations: `"median"` (median over iterations,
#'   the systematic response) or `"all"` (one line per iteration, the noise).
#' @return a `ggplot`.
#' @export
plot_length_indicators <- function(li, source = NULL, stat = c("median", "all")) {
  stat <- match.arg(stat)
  Lopt <- attr(li, "Lopt")
  if (!is.null(source)) li <- li[li$source %in% source, ]
  runs <- !is.null(li$run)
  if (stat == "median" && length(unique(li$iter)) > 1) {
    by <- intersect(c("run", "year", "source"), names(li))
    vars <- intersect(c("Lmean", "Pmat", "Pmega", "Z_lbar", "Z_true"), names(li))
    li <- stats::aggregate(li[, vars], li[, by, drop = FALSE], stats::median)
    li$iter <- "median"
  }
  id <- intersect(c("run", "iter", "year", "source"), names(li))
  panels <- c(Lmean = "Mean length (cm)", Pmat = "Pmat", Pmega = "Pmega",
              Z_lbar = "Z from mean length")
  d <- do.call(rbind, lapply(names(panels), function(v)
    data.frame(li[, id], what = panels[[v]], value = li[[v]])))
  d$what <- factor(d$what, levels = panels)
  d$grp <- interaction(d$iter, d$source, if (runs) d$run else "", drop = TRUE)
  colour_by <- if (runs) "run" else "source"
  cols <- if (runs) {
    run_colours(levels(li$run))
  } else {
    srcs <- setdiff(unique(li$source), "survey")
    fleet_colours(srcs)
  }
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value, group = .data$grp,
                                       colour = .data[[colour_by]])) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::geom_hline(data = data.frame(what = factor(panels[["Lmean"]], levels = panels),
                                          y = Lopt),
                        ggplot2::aes(yintercept = .data$y), linetype = 3)
  # true Z (answer key) as a dashed line in the Z panel
  if (!is.null(li$Z_true)) {
    dz <- data.frame(li[, id], what = factor(panels[["Z_lbar"]], levels = panels),
                     value = li$Z_true)
    dz$grp <- interaction(dz$iter, dz$source, if (runs) dz$run else "", drop = TRUE)
    g <- g + ggplot2::geom_line(data = dz, linetype = 2, linewidth = 0.6)
  }
  g <- g + ggplot2::scale_colour_manual(values = cols)
  g <- g + ggplot2::labs(x = "Year", y = NULL, colour = NULL,
                         subtitle = paste0("Dotted: optimum length Lopt. Dashed: true Z (Fbar + M)",
                                           if (identical(unique(li$iter), "median"))
                                             ". Median over iterations")) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  if (runs) g + ggplot2::facet_grid(what ~ source, scales = "free_y") else
    g + ggplot2::facet_wrap(~what, scales = "free_y")
}

#' Indicator dashboard
#'
#' All series of one species on one page, each scaled to its own mean:
#' abundance indices (CPUE by fleet, survey), total catch, and mean length
#' and Pmega by source. Independent indicators of the same stock should
#' broadly agree; where they diverge, one of them is probably not tracking
#' the stock (e.g. effort creep, a gear change, a cohort passing through).
#'
#' @inheritParams plot_availability
#' @param truth add true SSB (scaled to its mean) to the abundance panel.
#' @return a `ggplot`.
#' @export
plot_indicators <- function(data, species = NULL, iter = 1, truth = FALSE) {
  species <- species %||% unique(data$catch$species)[1]
  rel <- function(v) v / mean(v, na.rm = TRUE)
  idx <- rbind(pick(data$cpue, species, iter), pick(data$survey, species, iter))
  d <- data.frame(panel = "Abundance indices", year = idx$year, source = idx$source,
                  value = idx$index)
  ct <- stats::aggregate(catch_t ~ year, pick(transform(data$catch, source = data$catch$fleet),
                                              species, iter), sum)
  d <- rbind(d, data.frame(panel = "Total catch", year = ct$year, source = "catch",
                           value = rel(ct$catch_t)))
  li <- length_indicators(data, species, iter)
  for (s in unique(li$source)) {
    z <- li[li$source == s, ]
    d <- rbind(d, data.frame(panel = "Mean length", year = z$year, source = s,
                             value = rel(z$Lmean)),
               data.frame(panel = "Large fish (Pmega)", year = z$year, source = s,
                          value = rel(z$Pmega + 1e-6)))
  }
  if (truth && !is.null(data$truth)) {
    # scale true SSB over the years with index data, as the indices are
    tr <- data$truth[data$truth$species == species & data$truth$iter == iter[1] &
                       data$truth$year >= min(idx$year), ]
    d <- rbind(d, data.frame(panel = "Abundance indices", year = tr$year,
                             source = "true SSB", value = rel(tr$ssb)))
  }
  d$panel <- factor(d$panel, levels = c("Abundance indices", "Total catch",
                                        "Mean length", "Large fish (Pmega)"))
  cols <- c(src_colours(data), catch = "#444444", `true SSB` = "#8A8A8A")
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value, colour = .data$source)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3, colour = "grey50") +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$source == "true SSB"),
                       linewidth = 0.7, show.legend = FALSE) +
    ggplot2::geom_point(size = 1.1) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~panel, scales = "free_y") +
    ggplot2::labs(x = "Year", y = "Relative to series mean", colour = NULL,
                  title = species,
                  subtitle = "Do the indicators tell the same story?") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}
