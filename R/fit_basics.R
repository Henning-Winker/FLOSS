# Basic fits to observed data: growth, length-weight, maturity, catch curves,
# natural mortality, selectivity and CPUE assumptions. Each returns a
# `floss_fit` that plot_fit() draws, with the true values added on request.

#' Fit von Bertalanffy growth to observed length-at-age
#'
#' Nonlinear least squares, \eqn{L = L_\infty(1 - e^{-k(a + \delta - t_0)})},
#' on the aged fish of the biological samples. Fish are caught during the
#' year, so the observed (integer) age is shifted by `age_offset` (0.5 =
#' mid-year) to make `t0` comparable with the true value.
#'
#' @param data a `floss_data` or a [read_dataset()] list.
#' @param species species name.
#' @param iter iteration(s) (pooled).
#' @param source sources to use; default all.
#' @param years years to use; default all.
#' @param t0 fix t0 at this value (`NULL` = estimate).
#' @param age_offset fraction of a year added to observed ages.
#' @return a `floss_fit` with `pars` (`linf`, `k`, `t0`).
#' @export
fit_growth <- function(data, species, iter = 1, source = NULL, years = NULL,
                       t0 = NULL, age_offset = 0.5) {
  b <- pick(data$bio, species, iter, source, years)
  b <- b[!is.na(b$age), ]
  d <- data.frame(age = b$age + age_offset, len = b$len_mm / 10,
                  source = b$source)
  linf0 <- stats::quantile(d$len, 0.99) * 1.05
  if (is.null(t0)) {
    m <- stats::nls(len ~ linf * (1 - exp(-k * (age - t0))), d,
                    start = list(linf = linf0, k = 0.3, t0 = 0),
                    algorithm = "port",
                    lower = c(0.5 * linf0, 0.01, -5), upper = c(3 * linf0, 3, 2))
    pars <- stats::coef(m)
  } else {
    tz <- t0
    m <- stats::nls(len ~ linf * (1 - exp(-k * (age - tz))), d,
                    start = list(linf = linf0, k = 0.3), algorithm = "port",
                    lower = c(0.5 * linf0, 0.01), upper = c(3 * linf0, 3))
    pars <- c(stats::coef(m), t0 = t0)
  }
  ag <- seq(0, max(d$age) + 1, length.out = 100)
  vb <- function(p, a) p[["linf"]] * (1 - exp(-p[["k"]] * (a - p[["t0"]])))
  tru <- true_lh(data, species)
  fit_obj("growth", species, pars, d,
          curve = data.frame(x = ag, y = vb(pars, ag)),
          truth = if (!is.null(tru)) c(linf = tru$linf, k = tru$k, t0 = tru$t0),
          truth_curve = if (!is.null(tru)) data.frame(
            x = ag, y = vb(c(linf = tru$linf, k = tru$k, t0 = tru$t0), ag)),
          model = m, n = nrow(d))
}

#' Fit the length-weight relationship
#'
#' Linear regression of \eqn{\log W} on \eqn{\log L} (W in g, L in cm), with
#' the lognormal bias correction \eqn{a = e^{\alpha + \sigma^2/2}}.
#' @inheritParams fit_growth
#' @return a `floss_fit` with `pars` (`a`, `b`).
#' @export
fit_lw <- function(data, species, iter = 1, source = NULL, years = NULL) {
  b <- pick(data$bio, species, iter, source, years)
  d <- data.frame(len = b$len_mm / 10, wt = b$wt_g, source = b$source)
  d <- d[d$len > 0 & d$wt > 0, ]
  m <- stats::lm(log(wt) ~ log(len), d)
  s2 <- summary(m)$sigma^2
  pars <- c(a = exp(stats::coef(m)[[1]] + s2 / 2), b = stats::coef(m)[[2]])
  L <- seq(min(d$len), max(d$len), length.out = 100)
  tru <- true_lh(data, species)
  fit_obj("lw", species, pars, d,
          curve = data.frame(x = L, y = pars[["a"]] * L^pars[["b"]]),
          truth = if (!is.null(tru)) c(a = tru$a, b = tru$b),
          truth_curve = if (!is.null(tru)) data.frame(x = L, y = tru$a * L^tru$b),
          model = m, n = nrow(d))
}

#' Fit a logistic maturity ogive
#'
#' Binomial GLM (logit link) of maturity (0/1) on length:
#' \eqn{L_{50} = -\beta_0/\beta_1}, \eqn{L_{95} = (\log 19 - \beta_0)/\beta_1}.
#' @inheritParams fit_growth
#' @param bin_cm length bin for the plotted observed proportions.
#' @return a `floss_fit` with `pars` (`L50`, `L95`).
#' @export
fit_maturity <- function(data, species, iter = 1, source = NULL, years = NULL,
                         bin_cm = NULL) {
  b <- pick(data$bio, species, iter, source, years)
  d <- data.frame(len = b$len_mm / 10, mature = b$mature, source = b$source)
  m <- stats::glm(mature ~ len, stats::binomial(), d)
  cf <- stats::coef(m)
  pars <- c(L50 = -cf[[1]] / cf[[2]], L95 = (log(19) - cf[[1]]) / cf[[2]])
  bin <- bin_cm %||% (lfd_bin_mm(data, species) / 10)
  d$bin <- (floor(d$len / bin) + 0.5) * bin
  L <- seq(min(d$len), max(d$len), length.out = 100)
  ogive <- function(L, L50, L95) 1 / (1 + exp(-log(19) * (L - L50) / (L95 - L50)))
  tru <- true_lh(data, species)
  tr <- if (!is.null(tru)) {
    m95 <- data$settings$oem$bio$mat_L95 %||% 1.15
    c(L50 = tru$lm, L95 = tru$lm * m95)
  }
  fit_obj("maturity", species, pars, d,
          curve = data.frame(x = L, y = stats::plogis(cf[[1]] + cf[[2]] * L)),
          truth = tr,
          truth_curve = if (!is.null(tr)) data.frame(x = L, y = ogive(L, tr[["L50"]], tr[["L95"]])),
          model = m, n = nrow(d))
}

#' Age-based catch curve
#'
#' Regression of \eqn{\log n_a} on age over the descending limb of the
#' catch-at-age (pooled over the selected years); total mortality
#' \eqn{Z = -\text{slope}}. By default the fit starts one age after the
#' mode and drops the last (plus-group) age.
#'
#' @param x catch-at-age from [caa_obs()], or a `floss_data` (then
#'   `caa_obs()` is called with `species`, `source`, `iter`, `years`).
#' @inheritParams fit_growth
#' @param source fleet for `caa_obs()`.
#' @param ages ages to include in the regression (default: descending limb).
#' @param om optional `floss_om` of the species: adds the true mean Z over
#'   the fitted ages and years.
#' @return a `floss_fit` with `pars` (`Z`) and, when the true M is known, `F
#'   = Z - M`.
#' @export
catch_curve <- function(x, species = NULL, source = NULL, iter = 1,
                        years = NULL, ages = NULL, om = NULL) {
  data <- NULL
  if (inherits(x, "floss_data") || !is.null(x$bio)) {
    data <- x
    x <- caa_obs(data, species, source, iter = iter, years = years)
  }
  species <- x$species[1]
  n <- stats::aggregate(n ~ age, x, sum)
  n <- n[order(n$age), ]
  if (is.null(ages)) {
    pk <- n$age[which.max(n$n)]
    ages <- n$age[n$age > pk & n$age < max(n$age) & n$n > 0]
  }
  d <- data.frame(age = n$age, logn = log(pmax(n$n, 1e-12)), used = n$age %in% ages,
                  n = n$n)
  d <- d[d$n > 0, ]
  m <- stats::lm(logn ~ age, d[d$used, ])
  Z <- -stats::coef(m)[[2]]
  tru <- if (!is.null(data)) true_lh(data, species)
  pars <- c(Z = Z)
  if (!is.null(tru)) pars <- c(pars, M_true = tru$M, F = Z - tru$M)
  truth <- NULL
  if (!is.null(om)) {
    yrs <- intersect(ac(unique(x$year)), dimnames(om$stock)$year)
    aa <- intersect(ac(ages), dimnames(om$stock)$age)
    zz <- (harvest(om$stock) + m(om$stock))[aa, yrs]
    truth <- c(Z = mean(c(zz)))
  }
  fit_obj("catch_curve", species, pars, d,
          curve = data.frame(x = range(ages), y = stats::predict(m, data.frame(age = range(ages)))),
          truth = truth, model = m, n = sum(d$used),
          extra = list(source = x$source[1], years = range(x$year)))
}

#' Length-converted catch curve and selectivity
#'
#' Pauly's length-converted catch curve on the pooled length samples of one
#' source, using growth parameters (e.g. from [fit_growth()]): relative age
#' \eqn{t(L) = t_0 - \log(1 - L/L_\infty)/k}, time to grow through a bin
#' \eqn{\Delta t}, and \eqn{\log(n/\Delta t)} against \eqn{t}. Z is minus the
#' slope of the descending limb. Selectivity on the ascending limb is the
#' ratio of observed to expected numbers (the descending line extrapolated
#' back), fitted with a logistic curve.
#'
#' @inheritParams fit_growth
#' @param source one fleet (or `"survey"`).
#' @param linf,k,t0 growth parameters (cm, per year).
#' @param lmax_prop bins above `lmax_prop` x Linf are excluded.
#' @return a `floss_fit` with `pars` (`Z`, `L50`, `L95`).
#' @export
catch_curve_length <- function(data, species, source, linf, k, t0 = 0,
                               iter = 1, years = NULL, lmax_prop = 0.95) {
  l <- pick(data$lfd, species, iter, source, years)
  if (!nrow(l)) stop("no length samples for this selection")
  bw <- lfd_bin_mm(data, species) / 10
  n <- stats::aggregate(n ~ len_mm, l, sum)
  L1 <- n$len_mm / 10
  L2 <- L1 + bw
  ok <- L2 < lmax_prop * linf
  n <- n[ok, ]; L1 <- L1[ok]; L2 <- L2[ok]
  Lm <- (L1 + L2) / 2
  dt <- log((linf - L1) / (linf - L2)) / k
  tt <- t0 - log(1 - Lm / linf) / k
  y <- log(n$n / dt)
  pk <- which.max(y)
  desc <- seq_along(y) > pk & n$n > 0
  m <- stats::lm(y ~ tt, data.frame(y = y[desc], tt = tt[desc]))
  Z <- -stats::coef(m)[[2]]
  # selectivity on the ascending limb
  asc <- seq_len(pk)
  P <- pmin(1, (n$n / dt)[asc] / exp(stats::predict(m, data.frame(tt = tt[asc]))))
  sel_pars <- c(L50 = NA_real_, L95 = NA_real_)
  if (length(asc) >= 3) {
    fs <- try(stats::nls(P ~ 1 / (1 + exp(-log(19) * (L - L50) / d95)),
                         data.frame(P = P, L = Lm[asc]),
                         start = list(L50 = Lm[asc][which.min(abs(P - 0.5))],
                                      d95 = max(bw, diff(range(Lm[asc])) / 3)),
                         algorithm = "port", lower = c(0, 0.01)), silent = TRUE)
    if (!inherits(fs, "try-error")) {
      cf <- stats::coef(fs)
      sel_pars <- c(L50 = cf[["L50"]], L95 = cf[["L50"]] + cf[["d95"]])
    }
  }
  d <- data.frame(t = tt, y = y, used = desc, len = Lm, P = NA_real_)
  d$P[asc] <- P
  # true selectivity at length, if the species spec is known
  tsel <- NULL
  sp <- true_sp(data, species)
  if (!is.null(sp) && source %in% names(sp$sel)) {
    s <- sp$sel[[source]]
    tsel <- switch(s$type,
      logistic = FLRef::sel_logistic_len(Lm, L50 = s$L50, L95 = s$L95),
      normal = FLRef::sel_normal_len(Lm, Lpeak = s$Lpeak, sd = s$sd),
      dnormal = FLRef::sel_dnormal_len(Lm, Lpeak = s$Lpeak, sd_left = s$sd_left,
                                       sd_right = s$sd_right))
    tsel <- data.frame(x = Lm, y = tsel / max(tsel))
  }
  fit_obj("catch_curve_length", species, c(Z = Z, sel_pars), d,
          curve = data.frame(x = range(tt[desc]),
                             y = stats::predict(m, data.frame(tt = range(tt[desc])))),
          truth = if (!is.null(sp) && source %in% names(sp$sel) &&
                      sp$sel[[source]]$type == "logistic")
            c(L50 = sp$sel[[source]]$L50, L95 = sp$sel[[source]]$L95),
          truth_curve = tsel, model = m, n = sum(desc),
          extra = list(source = source, linf = linf, k = k, t0 = t0))
}

#' Empirical natural mortality estimators
#'
#' @param linf,k growth parameters (cm, per year).
#' @param tmax maximum observed age (years).
#' @param temp mean water temperature (deg C) for Pauly's (1980) estimator;
#'   omitted when `NULL`.
#' @return a data frame with `method`, `M` and `inputs`.
#' @export
#' @examples
#' M_empirical(linf = 71.4, k = 0.21, tmax = 16)
M_empirical <- function(linf = NULL, k = NULL, tmax = NULL, temp = NULL) {
  out <- list()
  if (!is.null(tmax)) {
    out[["Then et al. 2015 (tmax)"]] <- c(4.899 * tmax^-0.916, "tmax")
    out[["Hoenig 1983 (tmax)"]] <- c(exp(1.46 - 1.01 * log(tmax)), "tmax")
  }
  if (!is.null(linf) && !is.null(k))
    out[["Then et al. 2015 (growth)"]] <- c(4.118 * k^0.73 * linf^-0.33, "linf, k")
  if (!is.null(linf) && !is.null(k) && !is.null(temp))
    out[["Pauly 1980 (growth, temperature)"]] <- c(
      10^(-0.0066 - 0.279 * log10(linf) + 0.6543 * log10(k) + 0.4634 * log10(temp)),
      "linf, k, temp")
  data.frame(method = names(out), M = as.numeric(vapply(out, `[`, "", 1)),
             inputs = vapply(out, `[`, "", 2), row.names = NULL)
}

#' Check the CPUE assumptions against the truth (trainers)
#'
#' CPUE is assumed proportional to vulnerable biomass, \eqn{C/E = qB}. Fits
#' \eqn{\log I = \log q + \beta \log B} (\eqn{\beta = 1}: proportional;
#' \eqn{< 1}: hyperstable; \eqn{> 1}: hyperdepletion) and the trend in
#' \eqn{\log(I/B^\beta)} over years (catchability drift, e.g. effort creep).
#' Needs the answer key (`index_true`), so it works on a `floss_data`.
#'
#' @inheritParams fit_growth
#' @param source fleet.
#' @return a `floss_fit` with `pars` (`beta`, `q_trend_pct` per year, `r`).
#' @export
cpue_check <- function(data, species, source, iter = 1) {
  if (is.null(data$cpue$index_true))
    stop("cpue_check() needs the true index: use the floss_data, not an export")
  d <- pick(data$cpue, species, iter, source)
  d <- d[d$index > 0 & d$index_true > 0, ]
  m <- stats::lm(log(index) ~ log(index_true), d)
  beta <- stats::coef(m)[[2]]
  d$logq <- log(d$index) - beta * log(d$index_true)
  mq <- stats::lm(logq ~ year, d)
  pars <- c(beta = beta, q_trend_pct = 100 * (exp(stats::coef(mq)[[2]]) - 1),
            r = stats::cor(log(d$index), log(d$index_true)))
  x <- seq(min(log(d$index_true)), max(log(d$index_true)), length.out = 50)
  fit_obj("cpue", species, pars, d,
          curve = data.frame(x = exp(x), y = exp(stats::predict(m, data.frame(index_true = exp(x))))),
          truth = NULL, model = m, n = nrow(d),
          extra = list(source = source, qmodel = mq))
}

#' Plot a basic fit
#'
#' Draws the observations, the fitted curve and the estimates; with
#' `truth = TRUE` (trainers) the true curve or values are added as dashed
#' lines.
#' @param fit a `floss_fit`.
#' @param truth add the truth where known.
#' @param what for `catch_curve_length()`: `"curve"` or `"sel"`
#'   (selectivity); for `cpue_check()`: `"scatter"` or `"time"`.
#' @return a `ggplot`.
#' @export
plot_fit <- function(fit, truth = FALSE, what = NULL) {
  stopifnot(inherits(fit, "floss_fit"))
  d <- fit$data
  p <- fit$pars
  lab <- paste(names(p), signif(p, 3), sep = " = ", collapse = "   ")
  if (truth && !is.null(fit$truth))
    lab <- paste0(lab, "\ntrue: ", paste(names(fit$truth), signif(fit$truth, 3),
                                         sep = " = ", collapse = "   "))
  col <- "#0A5F78"
  th <- ggplot2::theme_bw()
  add_truth <- function(g) {
    if (truth && !is.null(fit$truth_curve))
      g + ggplot2::geom_line(data = fit$truth_curve, ggplot2::aes(.data$x, .data$y),
                             linetype = 2, colour = "#A85400", linewidth = 0.8)
    else g
  }
  g <- switch(fit$type,
    growth = ggplot2::ggplot(d, ggplot2::aes(.data$age, .data$len)) +
      ggplot2::geom_jitter(width = 0.1, height = 0, alpha = 0.3, size = 1,
                           colour = "grey35") +
      ggplot2::geom_line(data = fit$curve, ggplot2::aes(.data$x, .data$y),
                         colour = col, linewidth = 1) +
      ggplot2::labs(x = "Age (years, mid-year)", y = "Length (cm)"),
    lw = ggplot2::ggplot(d, ggplot2::aes(.data$len, .data$wt / 1000)) +
      ggplot2::geom_point(alpha = 0.25, size = 1, colour = "grey35") +
      ggplot2::geom_line(data = data.frame(x = fit$curve$x, y = fit$curve$y / 1000),
                         ggplot2::aes(.data$x, .data$y), colour = col, linewidth = 1) +
      ggplot2::labs(x = "Length (cm)", y = "Weight (kg)"),
    maturity = {
      pm <- stats::aggregate(mature ~ bin, d, function(x) c(p = mean(x), n = length(x)))
      pm <- do.call(data.frame, pm)
      ggplot2::ggplot(pm, ggplot2::aes(.data$bin, .data$mature.p)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$mature.n), colour = "grey35",
                            alpha = 0.6) +
        ggplot2::scale_size_area(max_size = 5, name = "Fish") +
        ggplot2::geom_line(data = fit$curve, ggplot2::aes(.data$x, .data$y),
                           colour = col, linewidth = 1) +
        ggplot2::labs(x = "Length (cm)", y = "Proportion mature")
    },
    catch_curve = ggplot2::ggplot(d, ggplot2::aes(.data$age, .data$logn)) +
      ggplot2::geom_point(ggplot2::aes(shape = .data$used), size = 2.2,
                          colour = "grey25") +
      ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
      ggplot2::geom_line(data = fit$curve, ggplot2::aes(.data$x, .data$y),
                         colour = col, linewidth = 1) +
      ggplot2::labs(x = "Age", y = "log(numbers)",
                    caption = "Filled points: ages used in the regression"),
    catch_curve_length = if (identical(what, "sel")) {
      ds <- d[!is.na(d$P), ]
      L50 <- p[["L50"]]; L95 <- p[["L95"]]
      Lx <- seq(min(d$len), max(d$len), length.out = 100)
      g0 <- ggplot2::ggplot(ds, ggplot2::aes(.data$len, .data$P)) +
        ggplot2::geom_point(size = 2, colour = "grey25") +
        ggplot2::labs(x = "Length (cm)", y = "Selectivity (observed / expected)")
      if (is.finite(L50)) g0 <- g0 + ggplot2::geom_line(
        data = data.frame(x = Lx, y = 1 / (1 + exp(-log(19) * (Lx - L50) / (L95 - L50)))),
        ggplot2::aes(.data$x, .data$y), colour = col, linewidth = 1)
      g0
    } else {
      ggplot2::ggplot(d, ggplot2::aes(.data$t, .data$y)) +
        ggplot2::geom_point(ggplot2::aes(shape = .data$used), size = 2.2,
                            colour = "grey25") +
        ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
        ggplot2::geom_line(data = fit$curve, ggplot2::aes(.data$x, .data$y),
                           colour = col, linewidth = 1) +
        ggplot2::labs(x = "Relative age t(L) (years)", y = "log(n / dt)",
                      caption = "Filled points: length bins used for Z")
    },
    cpue = if (identical(what, "time")) {
      dd <- rbind(data.frame(year = d$year, v = d$index, s = "Observed CPUE"),
                  data.frame(year = d$year, v = d$index_true, s = "True vulnerable biomass"))
      ggplot2::ggplot(dd, ggplot2::aes(.data$year, .data$v, linetype = .data$s)) +
        ggplot2::geom_line(colour = col) + ggplot2::geom_point(colour = col, size = 1.4) +
        ggplot2::labs(x = "Year", y = "Relative to mean", linetype = NULL) +
        ggplot2::theme(legend.position = "bottom")
    } else {
      ggplot2::ggplot(d, ggplot2::aes(.data$index_true, .data$index)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey40") +
        ggplot2::geom_point(colour = "grey25", size = 2) +
        ggplot2::geom_line(data = fit$curve, ggplot2::aes(.data$x, .data$y),
                           colour = col, linewidth = 1) +
        ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::labs(x = "True vulnerable biomass (relative, log)",
                      y = "CPUE (relative, log)",
                      caption = "Dotted: proportional (beta = 1)")
    })
  if (!(fit$type == "catch_curve_length" && !identical(what, "sel")) &&
      fit$type != "cpue") g <- add_truth(g)
  if (fit$type == "catch_curve_length" && identical(what, "sel") && truth &&
      !is.null(fit$truth_curve))
    g <- g + ggplot2::geom_line(data = fit$truth_curve, ggplot2::aes(.data$x, .data$y),
                                linetype = 2, colour = "#A85400", linewidth = 0.8)
  src <- fit$extra$source
  g + ggplot2::labs(title = paste0(fit$species, if (!is.null(src)) paste0(" - ", src)),
                    subtitle = lab) + th
}

#' @exportS3Method base::print
print.floss_fit <- function(x, ...) {
  cat("FLOSS fit:", x$type, "-", x$species, " (n =", x$n, ")\n")
  print(signif(x$pars, 4))
  if (!is.null(x$truth)) {
    cat("true:\n")
    print(signif(x$truth, 4))
  }
  invisible(x)
}

# ---- helpers ----------------------------------------------------------------

fit_obj <- function(type, species, pars, data, curve = NULL, truth = NULL,
                    truth_curve = NULL, model = NULL, n = NA, extra = NULL) {
  structure(list(type = type, species = species, pars = pars, data = data,
                 curve = curve, truth = truth, truth_curve = truth_curve,
                 model = model, n = n, extra = extra), class = "floss_fit")
}

true_sp <- function(data, species) {
  sps <- data$settings$species
  if (is.null(sps)) return(NULL)
  hit <- Filter(function(s) identical(s$species, species), sps)
  if (length(hit)) hit[[1]] else NULL
}

true_lh <- function(data, species) {
  sp <- true_sp(data, species)
  if (!is.null(sp)) return(sp$lh)
  lh <- data$lhpars
  if (is.null(lh)) return(NULL)
  r <- lh[lh$species == species, ]
  if (!nrow(r)) return(NULL)
  list(linf = r$linf_cm, k = r$k, t0 = r$t0, a = r$a, b = r$b, lm = r$lm_cm,
       M = r$M, tmax = r$tmax)
}
