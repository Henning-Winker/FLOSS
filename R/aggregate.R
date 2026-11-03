# Aggregation of the observations, as an analyst might do it (for teaching):
# catch and effort by fleet or summed over fleets (and species), nominal
# CPUE (catch / effort), and length samples pooled over fleets, as sampled
# or raised to the catch. Naive aggregation ignores that fleets differ in
# catchability, selectivity and sampling, which is the point of the
# exercise: compare the aggregate with the fleet-level data and the truth.

#' Catch, effort and nominal CPUE: by fleet or aggregated
#'
#' * `catch_table()`: observed catch (t) by fleet, summed over fleets
#'   (`"aggregate"`), or summed over fleets and species (`"combined"`).
#' * `effort_table()`: observed effort (fishing days) by fleet or summed over
#'   fleets.
#' * `nominal_cpue()`: catch / effort (kg per fishing day) at the chosen
#'   level, and as an index relative to its mean. `"aggregate"` divides the
#'   catch of all fleets by the effort of all fleets, ignoring that a day of
#'   trap fishing and a day of linefishing catch differently;
#'   `"combined"` also sums the catch of all species.
#'
#' Nominal CPUE uses all years with catch and effort, unlike the CPUE
#' index of [sample_data()] (which is a standardised index for the sampled
#' years). With changing fleet mix or effort creep, nominal CPUE drifts away
#' from the stock: compare it with `cpue_true` (true catch / true effort)
#' and the true biomass.
#'
#' @param data a `floss_data`.
#' @param species species names; default all (for `"combined"`, the species
#'   summed).
#' @param by `"fleet"`, `"aggregate"` or `"combined"` (catch, CPUE).
#' @return a data frame with `species`, `iter`, `year`, `source` (fleet or
#'   `"all fleets"`) and the values (observed and true).
#' @export
#' @examples
#' \donttest{
#' x <- simulate_fishery(nsim = 2, iters = 1)
#' head(nominal_cpue(x$data, by = "aggregate"))
#' plot_aggregate(x$data, "Lutjanus bohar", what = "cpue", by = "fleet")
#' }
catch_table <- function(data, species = NULL, by = c("fleet", "aggregate", "combined")) {
  by <- match.arg(by)
  d <- data$catch
  if (!is.null(species)) d <- d[d$species %in% species, ]
  d$source <- d$fleet
  if (by == "fleet") return(d[, c("species", "iter", "year", "source", "catch_t", "catch_true")])
  if (by == "combined") d$species <- "all species"
  out <- stats::aggregate(cbind(catch_t, catch_true) ~ species + iter + year, d, sum)
  out$source <- "all fleets"
  out[order(out$species, out$iter, out$year), c("species", "iter", "year", "source", "catch_t", "catch_true")]
}

#' @rdname catch_table
#' @export
effort_table <- function(data, by = c("fleet", "aggregate")) {
  by <- match.arg(by)
  d <- data$effort
  d$source <- d$fleet
  if (by == "fleet") return(d[, c("iter", "year", "source", "effort", "effort_true")])
  out <- stats::aggregate(cbind(effort, effort_true) ~ iter + year, d, sum)
  out$source <- "all fleets"
  out[order(out$iter, out$year), c("iter", "year", "source", "effort", "effort_true")]
}

#' @rdname catch_table
#' @export
nominal_cpue <- function(data, species = NULL, by = c("fleet", "aggregate", "combined")) {
  by <- match.arg(by)
  ct <- catch_table(data, species, by)
  ef <- effort_table(data, if (by == "fleet") "fleet" else "aggregate")
  d <- merge(ct, ef, by = c("iter", "year", "source"))
  d <- d[d$effort > 0, ]
  d$cpue <- 1000 * d$catch_t / d$effort                 # kg per fishing day
  d$cpue_true <- 1000 * d$catch_true / d$effort_true
  g <- interaction(d$species, d$iter, d$source, drop = TRUE)
  d$index <- d$cpue / stats::ave(d$cpue, g, FUN = mean)
  d$index_true <- d$cpue_true / stats::ave(d$cpue_true, g, FUN = mean)
  d <- d[order(d$species, d$source, d$iter, d$year), ]
  rownames(d) <- NULL
  d[, c("species", "iter", "year", "source", "catch_t", "effort", "cpue", "index",
        "catch_true", "effort_true", "cpue_true", "index_true")]
}

#' Length samples: by fleet or pooled, as sampled or raised to the catch
#'
#' * `raise = "sample"`: the counts as measured. Pooled over fleets
#'   (`by = "aggregate"`), each fleet weighs by its sample size, not by its
#'   catch: a fleet that is sampled heavily but catches little dominates the
#'   pooled distribution.
#' * `raise = "catch"`: each fleet's (season-weighted) sample is raised to the
#'   fleet's observed catch: numbers = sample numbers x catch (t) / sample
#'   weight (t), with weights from the length-weight relationship. Pooled
#'   over fleets, each fleet then weighs by its catch.
#'
#' Seasons are pooled for the year: counts summed for `"sample"`, equal
#' season weights for `"catch"` ([lfd_annual()]). The survey is not included.
#'
#' @inheritParams catch_table
#' @param by `"fleet"` or `"aggregate"` (all fleets pooled, source
#'   `"all fleets"`).
#' @param raise `"sample"` or `"catch"`.
#' @return a data frame like `data$lfd` without `season` (`n` is numeric when
#'   raised).
#' @export
lfd_aggregate <- function(data, species = NULL, by = c("fleet", "aggregate"),
                          raise = c("sample", "catch")) {
  by <- match.arg(by)
  raise <- match.arg(raise)
  l <- data$lfd[data$lfd$source != "survey", ]
  if (!is.null(species)) l <- l[l$species %in% species, ]
  l <- lfd_annual(l, if (raise == "catch") "equal" else "none")
  if (raise == "catch") {
    lh <- data$lhpars
    l$n <- as.numeric(l$n)
    parts <- lapply(split(l, l$species), function(d) {
      sp <- d$species[1]
      bw <- lfd_bin_mm(data, sp) / 10
      a <- lh$a[lh$species == sp]; b <- lh$b[lh$species == sp]
      L <- d$len_mm / 10 + bw / 2
      d$w_t <- d$n * a * L^b / 1e6                     # sample weight (t)
      key <- interaction(d$iter, d$year, d$source, drop = TRUE)
      sw <- stats::ave(d$w_t, key, FUN = sum)
      ct <- data$catch[data$catch$species == sp, ]
      cc <- ct$catch_t[match(paste(d$iter, d$year, d$source), paste(ct$iter, ct$year, ct$fleet))]
      d$n <- ifelse(sw > 0 & !is.na(cc), d$n * cc / sw, 0)
      d$w_t <- NULL
      d
    })
    l <- do.call(rbind, parts)
  }
  if (by == "aggregate") {
    l <- stats::aggregate(n ~ species + iter + year + len_mm, l, sum)
    l$source <- "all fleets"
  }
  l <- l[order(l$species, l$iter, l$year, l$source, l$len_mm),
         c("species", "iter", "year", "source", "len_mm", "n")]
  if (raise == "catch") l$n <- round(l$n, 1)
  rownames(l) <- NULL
  attr(l, "raise") <- raise
  l
}

#' Plot catch, effort or nominal CPUE by fleet or aggregated
#'
#' One panel per fleet (or one for all fleets); with several iterations the
#' lines are coloured by iteration. With `truth = TRUE`, true values are
#' dashed: for nominal CPUE, true catch / true effort, and in grey the true
#' spawning biomass relative to its mean, to show how far nominal CPUE is
#' from the stock.
#'
#' @inheritParams catch_table
#' @param species species name (ignored for `by = "combined"`).
#' @param iter iteration(s).
#' @param what `"catch"`, `"effort"` or `"cpue"` (nominal).
#' @param truth show true values.
#' @return a `ggplot`.
#' @export
plot_aggregate <- function(data, species = NULL, iter = 1, what = c("cpue", "catch", "effort"),
                           by = c("fleet", "aggregate", "combined"), truth = FALSE) {
  what <- match.arg(what)
  by <- match.arg(by)
  species <- species %||% unique(data$catch$species)[1]
  if (what == "effort" && by == "combined") by <- "aggregate"
  sp <- if (by == "combined") NULL else species
  d <- switch(what, catch = catch_table(data, sp, by), effort = effort_table(data, by),
              cpue = nominal_cpue(data, sp, by))
  vcol <- switch(what, catch = "catch_t", effort = "effort", cpue = "index")
  d$value <- d[[vcol]]
  d$true <- d[[paste0(sub("_t$", "", vcol), "_true")]]
  d <- d[d$iter %in% iter, ]
  lev <- c(intersect(names(src_colours(data)), unique(d$source)), setdiff(unique(d$source), names(src_colours(data))))
  d$source <- factor(d$source, levels = lev)
  d$iter <- factor(d$iter)
  multi <- nlevels(d$iter) > 1
  cols <- if (multi) iter_colours(levels(d$iter)) else src_colours(data)
  ylab <- switch(what, catch = "Catch (t)", effort = "Effort (fishing days)",
                 cpue = "Nominal CPUE (relative to mean)")
  ttl <- paste(if (by == "combined" && what != "effort") "All species" else if (what == "effort") "Effort" else species,
               "-", switch(by, fleet = "by fleet", aggregate = "all fleets summed",
                           combined = "all species and fleets summed"))
  g <- ggplot2::ggplot(d, ggplot2::aes(.data$year, .data$value, group = .data$iter,
                                       colour = if (multi) .data$iter else .data$source)) +
    ggplot2::geom_line(linewidth = 0.5) + ggplot2::geom_point(size = 1.1) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~source, scales = "free_y") +
    ggplot2::labs(x = "Year", y = ylab, title = ttl, colour = "Iteration") +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = if (multi) "bottom" else "none")
  if (truth) {
    g <- g + ggplot2::geom_line(ggplot2::aes(y = .data$true), linetype = 2, linewidth = 0.4)
    if (what == "cpue" && by != "combined" && !is.null(data$truth)) {
      tr <- data$truth[data$truth$species == species & data$truth$iter %in% iter, ]
      tr$value <- tr$ssb / stats::ave(tr$ssb, tr$iter, FUN = mean)
      tr <- merge(tr, data.frame(source = levels(d$source)), by = NULL)
      tr$source <- factor(tr$source, levels = levels(d$source))
      g <- g + ggplot2::geom_line(data = tr, ggplot2::aes(.data$year, .data$value,
                                                         group = .data$iter),
                                  inherit.aes = FALSE, colour = "grey55", linewidth = 0.6) +
        ggplot2::labs(caption = "Dashed: true catch / true effort. Grey: true spawning biomass (relative to mean)")
    }
  }
  g
}
