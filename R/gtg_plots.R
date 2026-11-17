#' True catch-at-length of a growth-morph operating model
#'
#' Catch numbers by length, year and season for each fleet, from the growth morphs in the catch at age
#' (the same keys the length sampler uses). This is the length composition the samples are drawn from,
#' without sampling error: use it to check samples, or as a true `FLStockLen` input for length-based methods.
#' @param om a `floss_om` built with `gtg`.
#' @param oem optional [oem_settings()], for the length bins (default `oem_settings(om$fishery)`).
#' @param fleets fleets (default all).
#' @param seasons number of seasons in the year (equal length; catch split by the Baranov fraction).
#' @param iters iterations.
#' @param years years (default all OM years).
#' @param season_weights `"catch"` (default: each season's share of the catch), `"equal"` (each season
#'   rescaled to the same total, as when the samples are split equally over the seasons) or a numeric vector
#'   of season weights. Totals over the year are kept.
#' @return `FLQuants` by fleet, each with dimensions length (lower bin limit, cm) x year x 1 x season x 1 x iter.
#' @export
gtg_lfq <- function(om, oem = NULL, fleets = NULL, seasons = 12, iters = 1, years = NULL,
                    season_weights = "catch") {
  if (is.null(om$gtg)) stop("om was built without growth morphs (gtg)")
  G <- om$gtg$G
  oem <- oem %||% oem_settings(om$fishery, survey = NULL)
  fleets <- fleets %||% om$fishery$fleets
  gears <- floss_gears(om, oem)
  stk <- om$stock
  years <- as.character(years %||% dimnames(stk)$year)
  yi <- match(years, as.character(om$gtg$years_full))
  cn <- catch_fleet(om, "n")
  tt <- season_mid(seasons)
  out <- lapply(fleets, function(f) {
    lens <- as.numeric(colnames(gears[[f]]$condALK))
    bin <- if (length(lens) > 1) diff(lens)[1] else 1
    arr <- array(0, c(length(lens), length(years), 1, seasons, 1, length(iters)))
    for (j in seq_along(iters)) {
      K <- gtg_keys(G, gtg_catch_morph(om, iters[j])[[f]], tt, lens, bin, f)   # age x len x year x season
      for (y in seq_along(years)) {
        ca <- c(cn[[f]][, years[y], , , , iters[j]])
        z <- c(harvest(stk)[, years[y], , , , iters[j]] + m(stk)[, years[y], , , , iters[j]])
        fr <- season_catch_frac(z, seasons)
        for (s in seq_len(seasons)) arr[, y, 1, s, 1, j] <- crossprod(K[, , yi[y], s], ca * fr[, s])
        if (!identical(season_weights, "catch")) {
          w <- if (identical(season_weights, "equal")) rep(1, seasons) else season_weights
          tot <- colSums(matrix(arr[, y, 1, , 1, j], ncol = seasons))
          sc <- ifelse(tot > 0, w / sum(w) * sum(tot) / tot, 0)
          arr[, y, 1, , 1, j] <- sweep(matrix(arr[, y, 1, , 1, j], ncol = seasons), 2, sc, `*`)
        }
      }
    }
    FLQuant(arr, dimnames = list(len = lens, year = years, unit = "unique", season = seq_len(seasons),
                                 area = "unique", iter = iters), units = "1000")
  })
  FLQuants(stats::setNames(out, fleets))
}

#' Plots of the growth-morph structure of an operating model
#'
#' * `plot_gtg_linf()`: distribution of asymptotic length. Each morph contributes a normal component
#'   (SD `cv_res` x its Linf) weighted by its recruitment share; the mixture matches the target normal
#'   with the total CV.
#' * `plot_gtg_growth()`: growth curves of the morphs (line width by recruitment share) and the mean VB curve.
#' * `plot_gtg_selectivity()`: effective selectivity at age by fleet (fleet F at age / apical F) in the
#'   first and last years, against the selectivity at age of the age-based OM (`FLRef::sel_la`).
#'   Removing fast growers first lowers effective selectivity at older ages as fishing goes on.
#' * `plot_gtg_lenage()`: mean length at age in the population and in each fleet's catch in one year,
#'   against the VB curve (Rosa Lee effect in the population; gear selection in the catch).
#' * `plot_gtg_lfq()`: true catch-at-length by fleet (annual, summed over seasons), with the sampled length
#'   compositions on top when `data` is given (both as proportions).
#' * `plot_gtg_spr()`: true equilibrium SPR by year from the morphs, with SPR at MSY of the morph equilibrium.
#' @param om a `floss_om` built with `gtg`.
#' @param iter iteration (`iters` for ribbons).
#' @param year year (default last OM year).
#' @param years years to show (default the last four).
#' @param data optional output of [sample_data()] or a data frame of length samples (`sample_lfd()`).
#' @param oem,seasons,season_weights passed to [gtg_lfq()]; use `season_weights = "equal"` to compare with
#'   samples split equally over the seasons.
#' @return a `ggplot`.
#' @name plot_gtg
NULL

#' @rdname plot_gtg
#' @export
plot_gtg_linf <- function(om) {
  G <- gtg_G(om); st <- G$settings; linf <- sum(G$recP * G$linf_g)
  x <- seq(linf * (1 - 4 * st$cv), linf * (1 + 4 * st$cv), length.out = 400)
  comp <- do.call(rbind, lapply(seq_along(G$linf_g), function(g)
    data.frame(x = x, morph = g, y = G$recP[g] * stats::dnorm(x, G$linf_g[g], st$cv_res * G$linf_g[g]))))
  mix <- stats::aggregate(y ~ x, comp, sum)
  tgt <- data.frame(x = x, y = stats::dnorm(x, linf, st$cv * linf))
  ggplot2::ggplot(comp, ggplot2::aes(x, y)) +
    ggplot2::geom_area(ggplot2::aes(group = morph, fill = factor(morph)), position = "identity",
                       alpha = 0.35, colour = NA) +
    ggplot2::geom_line(data = tgt, linewidth = 0.9, colour = "black") +
    ggplot2::geom_line(data = mix, linewidth = 0.9, colour = "grey55", linetype = 2) +
    ggplot2::scale_fill_manual(values = grDevices::hcl.colors(length(G$linf_g), "viridis"), guide = "none") +
    ggplot2::labs(x = expression(L[infinity] ~ "(cm)"), y = "Density",
                  subtitle = sprintf("%d morphs; CV total %.2f = between %.3f + within %.3f",
                                     length(G$linf_g), st$cv, st$cv_between, st$cv_res)) +
    ggplot2::theme_bw()
}

#' @rdname plot_gtg
#' @export
plot_gtg_growth <- function(om) {
  G <- gtg_G(om); lh <- om$species$lh
  t <- seq(0, max(G$ages) + 1, by = 0.1)
  L <- G$len(t)
  d <- data.frame(age = rep(t, each = nrow(L)), L = c(L), morph = rep(seq_len(nrow(L)), length(t)),
                  w = rep(G$recP, length(t)))
  vb <- data.frame(age = t, L = lh$linf * (1 - exp(-lh$k * (t - lh$t0))))
  ggplot2::ggplot(d, ggplot2::aes(age, L, group = morph)) +
    ggplot2::geom_line(ggplot2::aes(colour = morph, linewidth = w)) +
    ggplot2::geom_line(data = vb, ggplot2::aes(group = 1), colour = "black", linetype = 2) +
    ggplot2::scale_colour_viridis_c(guide = "none") +
    ggplot2::scale_linewidth(range = c(0.2, 1.4), guide = "none") +
    ggplot2::labs(x = "Age (years)", y = "Length (cm)", title = om$species$species) +
    ggplot2::theme_bw()
}

#' @rdname plot_gtg
#' @export
plot_gtg_selectivity <- function(om, iter = 1) {
  yrs <- dimnames(om$stock)$year
  yy <- yrs[c(1, length(yrs))]
  fl <- names(om$f_fleet)
  d <- do.call(rbind, lapply(fl, function(f) do.call(rbind, lapply(yy, function(y) {
    s <- c(om$f_fleet[[f]][, y, , , , iter]); s <- s / max(s)
    data.frame(fleet = f, age = an(dimnames(om$stock)$age), sel = s, what = paste("morphs", y))
  }))))
  a <- do.call(rbind, lapply(fl, function(f)
    data.frame(fleet = f, age = an(dimnames(om$stock)$age), sel = c(om$sel[[f]]$sel_a), what = "age-based OM")))
  d <- rbind(a, d)
  d$what <- factor(d$what, levels = unique(d$what))
  ggplot2::ggplot(d, ggplot2::aes(age, sel, colour = what, linetype = what)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~fleet) +
    ggplot2::scale_colour_manual(values = c("grey40", "#1B9E77", "#D95F02")) +
    ggplot2::scale_linetype_manual(values = c(2, 1, 1)) +
    ggplot2::labs(x = "Age (years)", y = "Selectivity at age (F / max F)", colour = NULL, linetype = NULL,
                  title = om$species$species) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

#' @rdname plot_gtg
#' @export
plot_gtg_lenage <- function(om, iter = 1, year = NULL) {
  G <- om$gtg$G; lh <- om$species$lh
  year <- as.character(year %||% utils::tail(dimnames(om$stock)$year, 1))
  y <- match(year, as.character(om$gtg$years_full))
  Fy <- sapply(G$fleets, function(f) c(om$gtg$fapic_full[[f]][, , , , , iter]))
  tr <- gtg_track(G, matrix(Fy, ncol = length(G$fleets), dimnames = list(NULL, G$fleets)))
  cn <- gtg_aggregate(G, matrix(Fy, ncol = length(G$fleets), dimnames = list(NULL, G$fleets)), tr)$cn
  Lm <- G$len(G$ages + 0.5)
  ml <- function(w) colSums(w * Lm) / pmax(colSums(w), 1e-300)
  d <- rbind(data.frame(age = G$ages, L = ml(tr$l[, , y]), what = "population"),
             do.call(rbind, lapply(G$fleets, function(f)
               data.frame(age = G$ages, L = ml(cn[[f]][, , y]), what = f))))
  wc <- c(rep(1, length(G$ages)), unlist(lapply(G$fleets, function(f) {
    w <- colSums(cn[[f]][, , y]); w / max(w)})))
  d <- d[wc > 1e-3, ]                                      # drop ages a fleet hardly catches
  vb <- data.frame(age = G$ages, L = lh$linf * (1 - exp(-lh$k * (G$ages + 0.5 - lh$t0))))
  cols <- c(population = "black", fleet_colours(G$fleets)[G$fleets])
  ggplot2::ggplot(d, ggplot2::aes(age, L, colour = what)) +
    ggplot2::geom_line(data = vb, ggplot2::aes(age, L), inherit.aes = FALSE, colour = "grey60",
                       linetype = 2) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(x = "Age (years, mid-year)", y = "Mean length (cm)", colour = NULL,
                  title = paste(om$species$species, year), subtitle = "dashed: VB curve of the mean morph") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
}

#' @rdname plot_gtg
#' @export
plot_gtg_lfq <- function(om, data = NULL, iter = 1, years = NULL, oem = NULL, seasons = 12,
                         season_weights = "catch") {
  yrs <- dimnames(om$stock)$year
  years <- as.character(years %||% utils::tail(yrs, 4))
  lq <- gtg_lfq(om, oem, seasons = seasons, iters = iter, years = years, season_weights = season_weights)
  d <- do.call(rbind, lapply(names(lq), function(f) {
    x <- as.data.frame(seasonSums(lq[[f]]))
    data.frame(fleet = f, year = as.character(x$year), L = an(as.character(x$len)), n = x$data)
  }))
  d$p <- d$n / stats::ave(d$n, d$fleet, d$year, FUN = sum)
  bin <- diff(sort(unique(d$L)))[1]
  p <- ggplot2::ggplot(d, ggplot2::aes(L + bin / 2, p)) +
    ggplot2::facet_grid(year ~ fleet)
  if (!is.null(data)) {
    l <- if (is.data.frame(data)) data else data$lfd
    l <- l[l$species == om$species$species & l$iter == iter & as.character(l$year) %in% years &
             l$source %in% names(lq), ]
    s <- stats::aggregate(n ~ year + source + len_mm, l, sum)
    s <- data.frame(fleet = s$source, year = as.character(s$year), L = s$len_mm / 10, n = s$n)
    s$p <- s$n / stats::ave(s$n, s$fleet, s$year, FUN = sum)
    p <- p + ggplot2::geom_col(data = s, ggplot2::aes(fill = fleet), width = bin, alpha = 0.6) +
      ggplot2::scale_fill_manual(values = fleet_colours(names(lq)), guide = "none")
  }
  p + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(x = "Length (cm)", y = "Proportion", title = om$species$species,
                  subtitle = if (is.null(data)) "true catch-at-length" else "line: truth; bars: samples") +
    ggplot2::theme_bw()
}

#' @rdname plot_gtg
#' @export
plot_gtg_spr <- function(om, iters = 1:3) {
  s <- as.data.frame(om$gtg$spr)
  q <- stats::aggregate(data ~ year, s, function(v) stats::quantile(v, c(0.05, 0.5, 0.95)))
  q <- data.frame(year = an(as.character(q$year)), q$data)
  names(q)[2:4] <- c("lo", "med", "hi")
  it <- s[s$iter %in% as.character(iters), ]
  it$year <- an(as.character(it$year))
  ggplot2::ggplot(q, ggplot2::aes(year, med)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80") +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_line(data = it, ggplot2::aes(year, data, group = iter), linewidth = 0.3, colour = "#D95F02") +
    ggplot2::geom_hline(yintercept = om$gtg$refpts[["SPRmsy"]], linetype = 2, colour = "#1B9E77") +
    ggplot2::geom_hline(yintercept = 0.4, linetype = 3) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Year", y = "SPR (morph equilibrium)", title = om$species$species,
                  subtitle = sprintf("dashed: SPR at MSY = %.2f; dotted: 0.4", om$gtg$refpts[["SPRmsy"]])) +
    ggplot2::theme_bw()
}

gtg_G <- function(om) {
  if (inherits(om, "floss_om")) {
    if (is.null(om$gtg)) stop("om was built without growth morphs (gtg)")
    return(om$gtg$G)
  }
  om
}
