# Plots of observed data (no model fits).

# fleet colours in fishery order (as in the app), then the survey
src_colours <- function(data) {
  fleet_colours(data$settings$fishery$fleets %||% unique(data$catch$fleet))
}

#' Plot which data exist for which years
#'
#' One row per data stream and source, one tile per year with data.
#' @param data a `floss_data` (or [read_dataset()] list).
#' @param species species name; default the first.
#' @param iter iteration(s); several iterations are overlaid (lines) or
#'   shown in rows (length samples).
#' @return a `ggplot`.
#' @export
plot_availability <- function(data, species = NULL, iter = 1) {
  species <- species %||% unique(data$catch$species)[1]
  get <- function(df, stream, srccol = "source") {
    if (is.null(df) || !nrow(df)) return(NULL)
    if (!is.null(df$species)) df <- df[df$species == species, ]
    if (!is.null(df$iter)) df <- df[df$iter %in% iter, ]
    if (!nrow(df)) return(NULL)
    u <- unique(df[, c("year", srccol)])
    data.frame(stream = stream, source = u[[srccol]], year = u$year)
  }
  d <- rbind(get(data$catch, "Catch", "fleet"), get(data$effort, "Effort", "fleet"),
             get(data$cpue, "CPUE"), get(data$survey, "Survey index"),
             get(data$lfd, "Length frequencies"), get(data$bio, "Biological samples"))
  d$stream <- factor(d$stream, levels = c("Catch", "Effort", "CPUE",
    "Survey index", "Length frequencies", "Biological samples"))
  d$source <- factor(d$source, levels = rev(names(src_colours(data))))
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$source, fill = .data$source)) +
    ggplot2::geom_tile(width = 0.85, height = 0.7) +
    ggplot2::facet_grid(stream ~ ., scales = "free_y", space = "free_y",
                        switch = "y") +
    ggplot2::scale_fill_manual(values = src_colours(data)) +
    ggplot2::labs(x = "Year", y = NULL, title = species) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
                   panel.grid.minor = ggplot2::element_blank())
}

#' Plot observed catch, CPUE and survey indices
#'
#' @inheritParams plot_availability
#' @param what `"catch"` or `"index"` (CPUE and survey).
#' @param truth add the true values (dashed; trainer view).
#' @return a `ggplot`.
#' @export
plot_obs <- function(data, species = NULL, iter = 1, what = c("index", "catch"),
                     truth = FALSE) {
  what <- match.arg(what)
  species <- species %||% unique(data$catch$species)[1]
  if (what == "catch") {
    cdat <- data$catch
    cdat$source <- cdat$fleet
    d <- pick(cdat, species, iter)
    d$value <- d$catch_t
    d$true <- d$catch_true
    ylab <- "Catch (t)"
  } else {
    d <- rbind(pick(data$cpue, species, iter), pick(data$survey, species, iter))
    d$value <- d$index
    d$true <- d$index_true
    ylab <- "Index (relative to mean)"
  }
  d$source <- factor(d$source, levels = intersect(names(src_colours(data)), d$source))
  d$iter <- factor(d$iter)
  multi <- nlevels(d$iter) > 1
  # gears are in separate panels, so with several iterations colour shows the
  # iteration; with one iteration it shows the gear
  col_by <- if (multi) "iter" else "source"
  cols <- if (multi) iter_colours(levels(d$iter)) else src_colours(data)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value,
                                       colour = .data[[col_by]],
                                       group = .data$iter)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = if (multi) 1.1 else 1.4) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~source, scales = "free_y") +
    ggplot2::labs(x = "Year", y = ylab, title = species, colour = "Iteration") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = if (multi) "bottom" else "none")
  if (truth && !is.null(d$true))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$true, colour = .data[[col_by]]),
                                linetype = 2, linewidth = 0.4, alpha = 0.8,
                                data = d, inherit.aes = TRUE) +
      ggplot2::labs(caption = "Dashed: true values")
  p
}

#' Plot length samples
#'
#' * `type = "bars"` (default): the classic length-frequency histogram, one
#'   panel per year (rows) and source (columns), see below.
#' * `type = "ridge"`: one length distribution per year, drawn as a ridge
#'   along the time axis (as FLRef's `plot_lfd_ridgeline()`), one panel per
#'   source; with `cohorts = TRUE`, von Bertalanffy lines follow each year
#'   class, so modes can be followed along them.
#' * `type = "bubble"`: one row per year, bubble area = proportion of the
#'   year's sample at each length. With `cohorts = TRUE`, von Bertalanffy
#'   growth lines follow each year class through the samples: real modes
#'   move along these lines; patterns that jump or vanish are sampling noise
#'   or selectivity, not dynamics.
#' * `type = "bars"`: the classic length-frequency histogram, one panel per
#'   year (rows) and source (columns), as proportions of each year's sample.
#'   With `cohorts = TRUE`, the mean length of each age class in that year is
#'   marked at the top, so modes can be matched to year classes.
#'
#' Seasonal samples are pooled for the year (counts summed).
#'
#' @inheritParams plot_availability
#' @param source fleets / survey to show; default all.
#' @param type `"bars"`, `"ridge"` or `"bubble"`.
#' @param years years to show; default the last six for `"bars"`, all for
#'   `"ridge"` and `"bubble"`.
#' @param cohorts add cohort growth lines (bubble) or age-class marks (bars).
#' @param season seasons to include (pooled); default all.
#' @param lfd optional length table to plot instead of `data$lfd`, e.g. from
#'   [lfd_aggregate()] (pooled over fleets or raised to the catch).
#' @return a `ggplot`.
#' @export
#' @examples
#' \donttest{
#' x <- simulate_fishery("lutjanus_bohar", nsim = 2, iters = 1)
#' plot_lfd(x$data, cohorts = TRUE)
#' plot_lfd(x$data, type = "bars", source = c("gillnet", "linefish"))
#' }
plot_lfd <- function(data, species = NULL, iter = 1, source = NULL,
                     type = c("bars", "ridge", "bubble"), years = NULL,
                     cohorts = FALSE, season = NULL, lfd = NULL) {
  type <- match.arg(type)
  species <- species %||% unique(data$lfd$species)[1]
  l <- lfd %||% data$lfd
  raised <- identical(attr(lfd, "raise"), "catch")
  if (!is.null(season) && !is.null(l$season)) l <- l[l$season %in% season, ]
  d <- pick(lfd_annual(l), species, iter[1], source)
  if (type == "bars" && is.null(years)) years <- utils::tail(sort(unique(d$year)), 6)
  if (!is.null(years)) d <- d[d$year %in% years, ]
  if (!nrow(d)) stop("no length samples for this selection")
  d$p <- d$n / stats::ave(d$n, d$year, d$source, d$iter, FUN = sum)
  d$source <- factor(d$source, levels = intersect(names(src_colours(data)), d$source))
  bw <- lfd_bin_mm(data, species) / 10
  lh <- if (cohorts) true_lh(data, species)
  if (type == "ridge") {
    d$L <- d$len_mm / 10 + bw / 2
    yrs <- sort(unique(d$year))
    step <- if (length(yrs) > 1) min(diff(yrs)) else 1
    # each year's distribution scaled to 90% of the space to the next year
    d$dens <- d$p / stats::ave(d$p, d$year, d$source, FUN = max) * 0.9 * step
    d$grp <- interaction(d$source, d$year)
    d <- d[order(d$source, d$year, d$L), ]
    g <- ggplot2::ggplot(d) +
      ggplot2::geom_ribbon(ggplot2::aes(x = .data$L, ymin = .data$year, ymax = .data$year + .data$dens,
                                        group = .data$grp, fill = .data$source), alpha = 0.75, colour = NA) +
      ggplot2::scale_fill_manual(values = src_colours(data), guide = "none") +
      ggplot2::coord_flip() + ggplot2::facet_wrap(~source, ncol = 1) +
      ggplot2::labs(x = "Length (cm)", y = "Year", title = species,
                    subtitle = paste("Iteration", iter[1], if (raised) "- raised to the catch" else "")) +
      ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
    if (!is.null(lh)) {
      yr <- range(d$year) + c(0, step)
      ages <- seq(0, lh$tmax + 1, by = 0.05)
      co <- do.call(rbind, lapply(seq(floor(yr[1]) - lh$tmax - 1, ceiling(yr[2])), function(cy)
        data.frame(t = cy + ages, L = lh$linf * (1 - exp(-lh$k * (ages - lh$t0))), cohort = cy)))
      co <- co[co$t >= yr[1] & co$t <= yr[2] & co$L >= min(d$L) & co$L <= max(d$L), ]
      g <- g + ggplot2::geom_line(data = co, ggplot2::aes(x = .data$L, y = .data$t, group = .data$cohort),
                                  colour = "grey45", linetype = 2, linewidth = 0.3, inherit.aes = FALSE)
    }
    return(g)
  }
  if (type == "bars") {
    d$L <- d$len_mm / 10 + bw / 2
    g <- ggplot2::ggplot(d, ggplot2::aes(.data$L, .data$p, fill = .data$source)) +
      ggplot2::geom_col(width = bw * 0.95) +
      ggplot2::scale_fill_manual(values = src_colours(data), guide = "none") +
      ggplot2::facet_grid(year ~ source, scales = "free_y") +
      ggplot2::labs(x = "Length (cm)", y = if (raised) "Proportion of the year's catch (raised)" else "Proportion of the year's sample",
                    title = species, subtitle = paste("Iteration", iter[1])) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
    if (!is.null(lh)) {
      ages <- 0:lh$tmax
      La <- lh$linf * (1 - exp(-lh$k * (ages + 0.5 - lh$t0)))
      ma <- data.frame(age = ages, L = La)
      # older ages crowd near Linf: keep ages at least two bins apart
      ma <- ma[c(TRUE, diff(ma$L) >= 2 * bw), ]
      ma <- ma[ma$L >= min(d$L) & ma$L <= max(d$L), ]
      g <- g + ggplot2::geom_vline(data = ma, ggplot2::aes(xintercept = .data$L),
                                   colour = "grey55", linewidth = 0.25, linetype = 3) +
        ggplot2::labs(caption = "Dotted: mean length of each age class in mid-year (young ages, left)")
    }
    return(g)
  }
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$len_mm / 10, .data$year)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$p, colour = .data$source),
                        alpha = 0.55) +
    ggplot2::scale_size_area(max_size = 5, guide = "none") +
    ggplot2::scale_colour_manual(values = src_colours(data),
                                 guide = "none") +
    ggplot2::facet_wrap(~source) +
    ggplot2::labs(x = "Length (cm)", y = "Year", title = species,
                  subtitle = "Bubble area = proportion of the year's sample") +
    ggplot2::theme_bw()
  if (!is.null(lh)) {
    yrs <- seq(min(d$year), max(d$year))
    lmax <- max(d$len_mm) / 10
    co <- do.call(rbind, lapply(seq(min(yrs) - lh$tmax, max(yrs)), function(b) {
      age <- yrs - b
      L <- lh$linf * (1 - exp(-lh$k * (age + 0.5 - lh$t0)))
      keep <- age >= 0 & L <= lmax
      if (sum(keep) < 2) return(NULL)
      data.frame(cohort = b, year = yrs[keep], L = L[keep])
    }))
    g <- g + ggplot2::geom_path(data = co, ggplot2::aes(.data$L, .data$year,
                                                       group = .data$cohort),
                                colour = "grey45", linewidth = 0.25, alpha = 0.7,
                                inherit.aes = FALSE)
  }
  g
}

#' Plot biological samples (observations only)
#'
#' * `"age_length"`: length against observed age.
#' * `"length_weight"`: weight against length.
#' * `"maturity"`: observed proportion mature by length bin (point size =
#'   number of fish).
#' * `"alk"`: observed age-length key (counts; bubble area).
#' * `"length"`: length composition of the biological samples.
#'
#' Points are always coloured by fleet (or survey); `by` sets the panels.
#' With `truth = TRUE` the true relationship of the operating model is drawn
#' over the observations: the von Bertalanffy growth curve (at the middle of
#' the sampling season, or mid-year), the length-weight curve and the
#' maturity ogive.
#'
#' @inheritParams plot_availability
#' @param what one of the options above.
#' @param source sources to include; default all.
#' @param by panels: `"fleet"` (one panel per fleet), `"aggregate"` (all
#'   fleets in one panel, as a naive analyst would pool them) or `"season"`
#'   (one panel per sampling season: fish grow through the year, so length
#'   at age rises from season to season).
#' @param truth draw the true growth, length-weight and maturity curves.
#' @return a `ggplot`.
#' @export
plot_bio <- function(data, species = NULL, iter = 1,
                     what = c("age_length", "length_weight", "maturity", "alk",
                              "length"),
                     source = NULL, by = c("fleet", "aggregate", "season"),
                     truth = FALSE) {
  what <- match.arg(what)
  by <- match.arg(by)
  species <- species %||% unique(data$bio$species)[1]
  b <- pick(data$bio, species, iter, source)
  if (by == "season" && is.null(b$season)) by <- "fleet"
  cols <- src_colours(data)
  b$source <- factor(b$source, levels = intersect(names(cols), unique(b$source)))
  nS <- if (is.null(b$season)) 1 else max(b$season, na.rm = TRUE)
  if (by == "season") b$panel <- factor(paste("Season", b$season), levels = paste("Season", seq_len(nS)))
  fac <- switch(by, fleet = ggplot2::facet_wrap(~source),
                season = ggplot2::facet_wrap(~panel), aggregate = NULL)
  th <- ggplot2::theme_bw()
  lh <- if (truth) true_lh(data, species)
  # true curves, one per panel (season panels use the season mid-point)
  panels <- switch(by, fleet = data.frame(source = factor(levels(b$source), levels = levels(b$source))),
                   season = data.frame(panel = factor(levels(b$panel), levels = levels(b$panel)),
                                       t = (seq_len(nS) - 0.5) / nS),
                   aggregate = data.frame(dummy = 1))
  vb <- function(age, t) lh$linf * (1 - exp(-lh$k * (age + t - lh$t0)))
  true_layer <- function(df, x, y) {
    df <- merge(panels, df, by = NULL)
    ggplot2::geom_line(data = df, ggplot2::aes(.data[[x]], .data[[y]]), inherit.aes = FALSE,
                       colour = "black", linewidth = 0.8, linetype = 2)
  }
  cap <- if (!is.null(lh)) "Dashed: true relationship (operating model)"
  g <- switch(what,
    age_length = {
      b <- b[!is.na(b$age), ]
      p <- ggplot2::ggplot(b, ggplot2::aes(.data$age, .data$len_mm / 10, colour = .data$source)) +
        ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1) +
        ggplot2::labs(x = "Observed age (years)", y = "Length (cm)")
      if (!is.null(lh)) {
        ages <- seq(0, max(b$age, lh$tmax), length.out = 200)
        tm <- if (by == "season") panels$t else 0.5
        cur <- do.call(rbind, lapply(seq_along(tm), function(i) {
          d <- data.frame(age = ages, L = vb(ages, tm[i]))
          if (by == "season") d$panel <- panels$panel[i]
          d
        }))
        if (by == "season") {
          cur$panel <- factor(cur$panel, levels = levels(b$panel))
          p <- p + ggplot2::geom_line(data = cur, ggplot2::aes(.data$age, .data$L), inherit.aes = FALSE,
                                      colour = "black", linewidth = 0.8, linetype = 2)
        } else p <- p + true_layer(cur, "age", "L")
      }
      p
    },
    length_weight = {
      p <- ggplot2::ggplot(b, ggplot2::aes(.data$len_mm / 10, .data$wt_g / 1000, colour = .data$source)) +
        ggplot2::geom_point(alpha = 0.35, size = 1) +
        ggplot2::labs(x = "Length (cm)", y = "Weight (kg)")
      if (!is.null(lh)) {
        L <- seq(min(b$len_mm) / 10, max(b$len_mm) / 10, length.out = 200)
        p <- p + true_layer(data.frame(L = L, W = lh$a * L^lh$b / 1000), "L", "W")
      }
      p
    },
    maturity = {
      bin <- lfd_bin_mm(data, species)
      b$lb <- (floor(b$len_mm / bin) + 0.5) * bin / 10
      grp <- c("lb", "source", if (by == "season") "panel")
      m <- stats::aggregate(b$mature, b[grp], function(x) c(p = mean(x), n = length(x)))
      m <- do.call(data.frame, m)
      names(m)[names(m) == "x.p"] <- "p"; names(m)[names(m) == "x.n"] <- "n"
      p <- ggplot2::ggplot(m, ggplot2::aes(.data$lb, .data$p, colour = .data$source)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$n), alpha = 0.7) +
        ggplot2::scale_size_area(max_size = 5, name = "Fish") +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = "Length (cm)", y = "Proportion mature")
      if (!is.null(lh)) {
        L95 <- lh$lm * (data$settings$oem$bio$mat_L95 %||% 1.15)
        L <- seq(min(m$lb), max(m$lb), length.out = 200)
        og <- data.frame(L = L, P = 1 / (1 + exp(-log(19) * (L - lh$lm) / (L95 - lh$lm))))
        p <- p + true_layer(og, "L", "P")
      }
      p
    },
    alk = {
      k <- alk_obs(data, species, iter, source)
      d <- as.data.frame(as.table(k), responseName = "n")
      d <- d[d$n > 0, ]
      d$len <- as.numeric(as.character(d$len_mm)) / 10
      d$age <- as.integer(as.character(d$age))
      fac <- NULL
      ggplot2::ggplot(d, ggplot2::aes(.data$age, .data$len)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$n), colour = "#0A5F78", alpha = 0.6) +
        ggplot2::scale_size_area(max_size = 6, name = "Fish aged") +
        ggplot2::labs(x = "Age (years)", y = "Length bin (cm, lower limit)",
                      subtitle = "Observed age-length key, all fleets")
    },
    length = ggplot2::ggplot(b, ggplot2::aes(.data$len_mm / 10, fill = .data$source)) +
      ggplot2::geom_histogram(bins = 40, colour = NA, alpha = 0.85, position = "stack") +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(x = "Length (cm)", y = "Fish", fill = NULL))
  if (what == "alk") fac <- NULL
  if (what == "length" && by == "fleet") fac <- ggplot2::facet_wrap(~source, scales = "free_y")
  g + (if (what != "length" && what != "alk") ggplot2::scale_colour_manual(values = cols)) + fac +
    ggplot2::labs(title = species, colour = NULL, caption = if (what %in% c("age_length", "length_weight", "maturity")) cap) + th +
    ggplot2::theme(legend.position = "bottom")
}

#' Plot catch-at-age from [caa_obs()]
#' @param caa result of [caa_obs()].
#' @param type `"bubble"` (year x age) or `"bar"` (age composition by year).
#' @return a `ggplot`.
#' @export
plot_caa <- function(caa, type = c("bubble", "bar")) {
  type <- match.arg(type)
  ttl <- paste(caa$species[1], "-", caa$source[1], "catch-at-age from length samples x ALK")
  if (type == "bubble") {
    ggplot2::ggplot(caa, ggplot2::aes(.data$year, .data$age)) +
      ggplot2::geom_point(ggplot2::aes(size = .data$n), colour = "#0A5F78",
                          alpha = 0.6) +
      ggplot2::scale_size_area(max_size = 7, name = "Numbers") +
      ggplot2::labs(x = "Year", y = "Age", title = ttl) + ggplot2::theme_bw()
  } else {
    ggplot2::ggplot(caa, ggplot2::aes(.data$age, .data$n)) +
      ggplot2::geom_col(fill = "#0A5F78") +
      ggplot2::facet_wrap(~year) +
      ggplot2::labs(x = "Age", y = "Numbers", title = ttl) + ggplot2::theme_bw()
  }
}

#' Seasonal length sampling: sample sizes and mean length by season
#'
#' Left: fish measured per season and year (how balanced the sampling is).
#' Right: mean length by season and year, with the pooled annual mean
#' (black) and the mean with equal season weights (dashed). Fish grow
#' through the year, so later seasons have larger fish; when most samples
#' come from one part of the year, the pooled mean follows that season, and
#' a shift in the sampling season looks like a change in the stock.
#'
#' @inheritParams plot_lfd
#' @return a `ggplot`.
#' @export
plot_seasons <- function(data, species = NULL, iter = 1, source = NULL) {
  species <- species %||% unique(data$lfd$species)[1]
  l <- pick(data$lfd, species, iter[1], source)
  if (is.null(l$season)) stop("these length samples have no seasons")
  if (is.null(source)) l <- l[l$source != "survey", ]
  bw <- lfd_bin_mm(data, species) / 10
  l$L <- l$len_mm / 10 + bw / 2
  n <- stats::aggregate(n ~ year + season + source, l, sum)
  ml <- do.call(rbind, lapply(split(l, list(l$year, l$season, l$source), drop = TRUE), function(d)
    data.frame(year = d$year[1], season = d$season[1], source = d$source[1],
               value = stats::weighted.mean(d$L, d$n))))
  pool <- function(w) {
    a <- lfd_annual(l[, c("species", "iter", "year", "season", "source", "len_mm", "n")], w)
    a$L <- a$len_mm / 10 + bw / 2
    do.call(rbind, lapply(split(a, list(a$year, a$source), drop = TRUE), function(d)
      data.frame(year = d$year[1], source = d$source[1],
                 value = stats::weighted.mean(d$L, d$n), how = w)))
  }
  pm <- rbind(pool("none"), pool("equal"))
  pm$how <- factor(pm$how, levels = c("none", "equal"),
                   labels = c("pooled: counts summed", "pooled: equal season weights"))
  lev <- c("Fish measured", "Mean length (cm)")
  n$panel <- factor(lev[1], lev); ml$panel <- factor(lev[2], lev); pm$panel <- factor(lev[2], lev)
  n$value <- n$n
  n$season <- factor(n$season); ml$season <- factor(ml$season)
  lev_src <- intersect(names(src_colours(data)), unique(l$source))   # fishery order
  n$source <- factor(n$source, lev_src); ml$source <- factor(ml$source, lev_src)
  pm$source <- factor(pm$source, lev_src)
  ggplot2::ggplot(mapping = ggplot2::aes(.data$year, .data$value)) +
    ggplot2::geom_col(data = n, ggplot2::aes(fill = .data$season), width = 0.8) +
    ggplot2::geom_point(data = ml, ggplot2::aes(colour = .data$season), size = 1.7) +
    ggplot2::geom_line(data = pm, ggplot2::aes(linetype = .data$how), colour = "black") +
    ggplot2::facet_grid(panel ~ source, scales = "free_y", switch = "y") +
    ggplot2::scale_fill_manual(values = season_colours(nlevels(n$season)), aesthetics = c("colour", "fill")) +
    ggplot2::scale_linetype_manual(values = c(1, 2)) +
    ggplot2::labs(x = "Year", y = NULL, fill = "Season", colour = "Season", linetype = NULL,
                  title = species, subtitle = paste("Iteration", iter[1])) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                   strip.placement = "outside")
}
