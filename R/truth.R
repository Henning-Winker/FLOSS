# True stock status: the answer key summarised.

#' True stock status
#'
#' Summarises the true status of each species in each dataset (iteration):
#' B/Bmsy, F/Fmsy and catch/MSY in a given year, and the Kobe status
#' (`"Healthy"`: B >= Bmsy and F <= Fmsy; `"Overfishing"`: F > Fmsy with
#' B >= Bmsy; `"Overfished"`: B < Bmsy with F <= Fmsy; `"Overfished and
#' overfishing"`: both).
#'
#' @param x a `floss_data` (uses its answer key) or a `floss_om` /
#'   `floss_oms`.
#' @param year assessment year; default the last year.
#' @param iters iterations (for operating models; default all).
#' @param summary return the median over iterations per species instead.
#' @return a data frame.
#' @export
status_table <- function(x, year = NULL, iters = NULL, summary = FALSE) {
  tr <- truth_long(x, iters)
  year <- year %||% max(tr$year)
  d <- tr[tr$year == year, c("species", "iter", "year", "B_Bmsy", "F_Fmsy",
                             "C_MSY")]
  d$status <- kobe_status(d$B_Bmsy, d$F_Fmsy)
  d <- d[order(d$species, d$iter), ]
  rownames(d) <- NULL
  if (!summary) return(d)
  s <- stats::aggregate(cbind(B_Bmsy, F_Fmsy, C_MSY) ~ species + year, d,
                        stats::median)
  s$status <- kobe_status(s$B_Bmsy, s$F_Fmsy)
  s$p_healthy <- as.numeric(tapply(d$status == "Healthy", d$species, mean)[s$species])
  s
}

#' Plot true status trajectories
#'
#' B/Bmsy, F/Fmsy and catch/MSY over time by species, one line per dataset
#' (iteration), with the reference level 1 dashed (as `FLRef::plotAdvice()`
#' on ratios).
#' @inheritParams status_table
#' @param species species to show; default all.
#' @return a `ggplot`.
#' @export
plot_truth <- function(x, species = NULL, iters = NULL) {
  tr <- truth_long(x, iters)
  if (!is.null(species)) tr <- tr[tr$species %in% species, ]
  d <- rbind(
    data.frame(tr[, c("species", "iter", "year")], qname = "B/Bmsy", data = tr$B_Bmsy),
    data.frame(tr[, c("species", "iter", "year")], qname = "F/Fmsy", data = tr$F_Fmsy),
    data.frame(tr[, c("species", "iter", "year")], qname = "Catch/MSY", data = tr$C_MSY))
  d$qname <- factor(d$qname, levels = c("B/Bmsy", "F/Fmsy", "Catch/MSY"))
  ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$data,
                                  group = .data$iter)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
    ggplot2::geom_line(colour = "#0A5F78", alpha = 0.6, linewidth = 0.4) +
    ggplot2::facet_grid(qname ~ species, scales = "free_y") +
    ggplot2::labs(x = "Year", y = NULL,
                  subtitle = "True status, one line per dataset (iteration)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text.x = ggplot2::element_text(face = "italic"))
}

#' Kobe plot of true status
#'
#' True B/Bmsy against F/Fmsy, coloured by species. With several iterations,
#' one point per species and iteration in the assessment year. With one
#' iteration (or `trajectory = TRUE`), the whole history of each species is
#' drawn as a path from the first year (open circle) to the assessment year
#' (filled circle), as in JABBA's Kobe plot.
#'
#' @inheritParams status_table
#' @param species species to show; default all.
#' @param trajectory draw the history as a path; default when one iteration
#'   is shown.
#' @return a `ggplot`.
#' @export
plot_kobe <- function(x, year = NULL, iters = NULL, species = NULL, trajectory = NULL) {
  tr <- truth_long(x, iters)
  if (!is.null(species)) tr <- tr[tr$species %in% species, ]
  year <- year %||% max(tr$year)
  its <- unique(tr$iter)
  trajectory <- trajectory %||% (length(its) == 1)
  end <- tr[tr$year == year, ]
  path <- tr[tr$year <= year, ]
  path <- path[order(path$species, path$iter, path$year), ]
  shown <- if (trajectory) path else end
  mx <- max(2, shown$B_Bmsy, shown$F_Fmsy, na.rm = TRUE) * 1.05
  q <- data.frame(xmin = c(0, 1, 0, 1), xmax = c(1, mx, 1, mx),
                  ymin = c(0, 0, 1, 1), ymax = c(1, 1, mx, mx),
                  fill = c("#FFE08A", "#5CB85C", "#D9534F", "#F0A04B"))
  sps <- unique(tr$species)
  cols <- stats::setNames(floss_cols(length(sps)), sps)
  g <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = q, ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                              ymin = .data$ymin, ymax = .data$ymax),
                       fill = q$fill, alpha = 0.75)
  if (trajectory) {
    start <- path[!duplicated(paste(path$species, path$iter)), ]
    g <- g +
      ggplot2::geom_path(data = path, ggplot2::aes(.data$B_Bmsy, .data$F_Fmsy, colour = .data$species,
                                                   group = interaction(.data$species, .data$iter)),
                         linewidth = 0.7) +
      ggplot2::geom_point(data = start, ggplot2::aes(.data$B_Bmsy, .data$F_Fmsy, colour = .data$species),
                          shape = 21, fill = "white", size = 2.5, stroke = 1) +
      ggplot2::geom_point(data = end, ggplot2::aes(.data$B_Bmsy, .data$F_Fmsy, fill = .data$species),
                          shape = 21, colour = "black", size = 3.6, stroke = 0.8) +
      ggplot2::scale_fill_manual(values = cols, guide = "none")
    sub <- paste0("True status ", min(path$year), "-", year, ", iteration ",
                  paste(its, collapse = ", "), " (open: first year, filled: ", year, ")")
  } else {
    g <- g + ggplot2::geom_point(data = end, ggplot2::aes(.data$B_Bmsy, .data$F_Fmsy, colour = .data$species),
                                 size = 2.6, alpha = 0.9)
    sub <- paste0("True status in ", year, ", one point per species and iteration (",
                  length(its), " iterations)")
  }
  g + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_cartesian(xlim = c(0, mx), ylim = c(0, mx), expand = FALSE) +
    ggplot2::labs(x = "B / Bmsy", y = "F / Fmsy", colour = NULL, subtitle = sub) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.text = ggplot2::element_text(face = "italic"), legend.position = "bottom")
}

# ---- helpers ----------------------------------------------------------------

truth_long <- function(x, iters = NULL) {
  if (inherits(x, "floss_data")) {
    tr <- x$truth
    rp <- x$refpts
    if (!is.null(iters)) tr <- tr[tr$iter %in% iters, ]
  } else {
    oms <- if (inherits(x, "floss_om")) list(x) else x
    tr <- do.call(rbind, lapply(oms, function(o)
      truth_table(o, iters %||% seq_len(dims(o$stock)$iter))))
    rp <- do.call(rbind, lapply(oms, refpts_table))
  }
  tr$C_MSY <- tr$catch / rp$MSY[match(tr$species, rp$species)]
  tr
}

kobe_status <- function(b, f) {
  ifelse(b >= 1 & f <= 1, "Healthy",
         ifelse(b >= 1, "Overfishing",
                ifelse(f <= 1, "Overfished", "Overfished and overfishing")))
}

#' Download the truth (trainers)
#'
#' Writes the true time series (SSB, recruitment, Fbar, catch, B/Bmsy,
#' F/Fmsy, catch/MSY), the reference points and the status table, without
#' any observed data.
#' @inheritParams status_table
#' @param dir output folder (created).
#' @param format `"csv"` and/or `"xlsx"` (one workbook `truth.xlsx` with three
#'   sheets; needs 'writexl').
#' @return the folder, invisibly.
#' @export
export_truth <- function(x, dir, format = "csv", year = NULL, iters = NULL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tabs <- list(
    truth_timeseries = truth_long(x, iters),
    refpts = if (inherits(x, "floss_data")) x$refpts else
      do.call(rbind, lapply(if (inherits(x, "floss_om")) list(x) else x,
                            refpts_table)),
    status = status_table(x, year, iters))
  if ("csv" %in% format)
    for (n in names(tabs))
      utils::write.csv(tabs[[n]], file.path(dir, paste0(n, ".csv")),
                       row.names = FALSE)
  if ("xlsx" %in% format) {
    if (!requireNamespace("writexl", quietly = TRUE))
      stop("format = 'xlsx' needs the 'writexl' package")
    writexl::write_xlsx(tabs, file.path(dir, "truth.xlsx"))
  }
  invisible(dir)
}
