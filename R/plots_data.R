# Plots of observed data (no model fits).

# fleet colours in fishery order (as in the app); survey in navy
src_colours <- function(data) {
  fl <- data$settings$fishery$fleets %||% unique(data$catch$fleet)
  base <- c("#0A6E8C", "#A85400", "#5B45A8", "#2E7D32", "#8E2C48", "#555555")
  base <- rep(base, length.out = length(fl))
  c(stats::setNames(base, fl), survey = "#102330")
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
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value,
                                       colour = .data$source,
                                       group = .data$iter)) +
    ggplot2::geom_line(linewidth = 0.4, alpha = if (multi) 0.6 else 1) +
    ggplot2::geom_point(size = if (multi) 0.9 else 1.4,
                        alpha = if (multi) 0.6 else 1) +
    ggplot2::scale_colour_manual(values = src_colours(data)) +
    ggplot2::facet_wrap(~source, scales = "free_y") +
    ggplot2::labs(x = "Year", y = ylab, title = species,
                  subtitle = if (multi) paste("Iterations", paste(levels(d$iter),
                                                                   collapse = ", "))) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
  if (truth && !is.null(d$true))
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$true), linetype = 2,
                                colour = "grey20")
  p
}

#' Plot length samples as bubbles (year x length)
#'
#' @inheritParams plot_availability
#' @param source fleets / survey to show (facets); default all.
#' @return a `ggplot`.
#' @export
plot_lfd <- function(data, species = NULL, iter = 1, source = NULL) {
  species <- species %||% unique(data$lfd$species)[1]
  d <- pick(data$lfd, species, iter, source)
  d$p <- d$n / stats::ave(d$n, d$year, d$source, d$iter, FUN = sum)
  d$source <- factor(d$source, levels = intersect(names(src_colours(data)), d$source))
  ggplot2::ggplot(d, ggplot2::aes(.data$len_mm / 10, .data$year)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$p, colour = .data$source),
                        alpha = 0.55) +
    ggplot2::scale_size_area(max_size = 5, guide = "none") +
    ggplot2::scale_colour_manual(values = src_colours(data),
                                 guide = "none") +
    (if (length(unique(d$iter)) > 1)
      ggplot2::facet_grid(paste("Iteration", iter) ~ source) else
        ggplot2::facet_wrap(~source)) +
    ggplot2::labs(x = "Length (cm)", y = "Year", title = species,
                  subtitle = "Bubble area = proportion of the year's sample") +
    ggplot2::theme_bw()
}

#' Plot biological samples (observations only)
#'
#' * `"age_length"`: length against observed age.
#' * `"length_weight"`: weight against length.
#' * `"maturity"`: observed proportion mature by length bin (point size =
#'   number of fish).
#' * `"alk"`: observed age-length key (counts; bubble area).
#' * `"length"`: length composition of the biological samples by source.
#'
#' @inheritParams plot_availability
#' @param what one of the options above.
#' @param source sources to include; default all.
#' @return a `ggplot`.
#' @export
plot_bio <- function(data, species = NULL, iter = 1,
                     what = c("age_length", "length_weight", "maturity", "alk",
                              "length"),
                     source = NULL) {
  what <- match.arg(what)
  species <- species %||% unique(data$bio$species)[1]
  b <- pick(data$bio, species, iter, source)
  cols <- src_colours(data)
  th <- ggplot2::theme_bw()
  switch(what,
    age_length = {
      b <- b[!is.na(b$age), ]
      ggplot2::ggplot(b, ggplot2::aes(.data$age, .data$len_mm / 10,
                                      colour = .data$source)) +
        ggplot2::geom_jitter(width = 0.15, height = 0, alpha = 0.4, size = 1) +
        ggplot2::scale_colour_manual(values = cols) +
        ggplot2::labs(x = "Observed age (years)", y = "Length (cm)",
                      colour = NULL, title = species) + th
    },
    length_weight = ggplot2::ggplot(b, ggplot2::aes(.data$len_mm / 10,
                                                    .data$wt_g / 1000,
                                                    colour = .data$source)) +
      ggplot2::geom_point(alpha = 0.35, size = 1) +
      ggplot2::scale_colour_manual(values = cols) +
      ggplot2::labs(x = "Length (cm)", y = "Weight (kg)", colour = NULL,
                    title = species) + th,
    maturity = {
      bin <- lfd_bin_mm(data, species)
      b$lb <- (floor(b$len_mm / bin) + 0.5) * bin / 10
      m <- stats::aggregate(mature ~ lb, b, function(x) c(p = mean(x), n = length(x)))
      m <- do.call(data.frame, m)
      ggplot2::ggplot(m, ggplot2::aes(.data$lb, .data$mature.p)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$mature.n),
                            colour = "#0A5F78", alpha = 0.7) +
        ggplot2::scale_size_area(max_size = 5, name = "Fish") +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = "Length (cm)", y = "Proportion mature",
                      title = species) + th
    },
    alk = {
      k <- alk_obs(data, species, iter, source)
      d <- as.data.frame(as.table(k), responseName = "n")
      d <- d[d$n > 0, ]
      d$len <- as.numeric(as.character(d$len_mm)) / 10
      d$age <- as.integer(as.character(d$age))
      ggplot2::ggplot(d, ggplot2::aes(.data$age, .data$len)) +
        ggplot2::geom_point(ggplot2::aes(size = .data$n), colour = "#0A5F78",
                            alpha = 0.6) +
        ggplot2::scale_size_area(max_size = 6, name = "Fish aged") +
        ggplot2::labs(x = "Age (years)", y = "Length bin (cm, lower limit)",
                      title = paste(species, "- observed age-length key")) + th
    },
    length = ggplot2::ggplot(b, ggplot2::aes(.data$len_mm / 10,
                                             fill = .data$source)) +
      ggplot2::geom_histogram(bins = 40, colour = NA, alpha = 0.85) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::facet_wrap(~source, scales = "free_y") +
      ggplot2::labs(x = "Length (cm)", y = "Fish", title = species) + th +
      ggplot2::theme(legend.position = "none"))
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
