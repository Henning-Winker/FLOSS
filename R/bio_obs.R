# Observed biology from biological samples: age-length key and catch-at-age.
# Nothing is fitted: everything is counted or averaged from the observations.

#' Observed age-length key
#'
#' Counts aged fish by length bin and age, from the biological samples.
#'
#' @param data a `floss_data` (or a list with a `bio` and `lfd` data frame,
#'   e.g. from [read_dataset()]).
#' @param species species name.
#' @param iter iteration(s); several are pooled.
#' @param source sources (fleets / survey) to pool; default all.
#' @param years years to pool; default all.
#' @param bin_mm length bin (mm); default the bin of the length samples.
#' @param prop return \eqn{P(a \mid l)} (rows sum to one) instead of counts.
#' @return a matrix, length bins (lower limit, mm) x ages.
#' @export
alk_obs <- function(data, species, iter = 1, source = NULL, years = NULL,
                    bin_mm = NULL, prop = FALSE) {
  b <- pick(data$bio, species, iter, source, years)
  b <- b[!is.na(b$age), ]
  if (!nrow(b)) stop("no aged fish for this selection")
  if (is.null(bin_mm)) bin_mm <- lfd_bin_mm(data, species)
  lb <- floor(b$len_mm / bin_mm) * bin_mm
  tab <- table(len_mm = factor(lb, levels = seq(min(lb), max(lb), by = bin_mm)),
               age = factor(b$age, levels = seq(0, max(b$age))))
  m <- unclass(tab)
  if (prop) {
    rs <- rowSums(m)
    m <- m / ifelse(rs > 0, rs, NA)
  }
  m
}

#' Catch-at-age from length samples and an age-length key
#'
#' Applies an observed age-length key to the length samples of each year,
#' \eqn{\hat n_{a,y} = \sum_l n_{l,y} P(a \mid l)}. Length bins without aged
#' fish cannot be assigned an age and are dropped (their share is reported).
#' With `raise = "catch"` the numbers are raised so that their weight
#' (\eqn{W = a L^b}, g) equals the observed catch of the source (t).
#'
#' @inheritParams alk_obs
#' @param source one fleet (or `"survey"`).
#' @param alk `"pooled"` (one key over all years of the selection) or
#'   `"annual"` (a key per year).
#' @param raise `"none"` (sample numbers) or `"catch"`.
#' @return a data frame with `species`, `iter`, `year`, `source`, `age`,
#'   `n`, and attribute `unassigned` (share of measured fish per year in
#'   length bins without aged fish).
#' @export
caa_obs <- function(data, species, source, iter = 1, years = NULL,
                    alk = c("pooled", "annual"), raise = c("none", "catch")) {
  alk <- match.arg(alk)
  raise <- match.arg(raise)
  if (length(iter) != 1) stop("caa_obs() works on one iteration at a time")
  bin <- lfd_bin_mm(data, species)
  l <- pick(lfd_annual(data$lfd), species, iter, source, years)
  if (!nrow(l)) stop("no length samples for this selection")
  yrs <- sort(unique(l$year))
  key_pooled <- if (alk == "pooled")
    alk_obs(data, species, iter, source, yrs, bin, prop = TRUE)
  lh <- data$lhpars[data$lhpars$species == species, ]
  unassigned <- numeric()
  out <- lapply(yrs, function(y) {
    key <- if (alk == "pooled") key_pooled else
      tryCatch(alk_obs(data, species, iter, source, y, bin, prop = TRUE),
               error = function(e) NULL)
    if (is.null(key)) return(NULL)
    ly <- l[l$year == y, ]
    k <- key[match(ly$len_mm, as.numeric(rownames(key))), , drop = FALSE]
    ok <- stats::complete.cases(k)
    unassigned[as.character(y)] <<- sum(ly$n[!ok]) / sum(ly$n)
    n_a <- colSums(k[ok, , drop = FALSE] * ly$n[ok])
    if (raise == "catch") {
      cy <- data$catch[data$catch$species == species & data$catch$iter == iter &
                         data$catch$fleet == source & data$catch$year == y, ]
      w_l <- lh$a * ((ly$len_mm[ok] + bin / 2) / 10)^lh$b         # g
      sample_wt <- sum(ly$n[ok] * w_l)
      if (nrow(cy) && sample_wt > 0) n_a <- n_a * cy$catch_t * 1e6 / sample_wt
    }
    data.frame(species = species, iter = iter, year = y, source = source,
               age = as.integer(names(n_a)), n = unname(n_a))
  })
  out <- do.call(rbind, out)
  attr(out, "unassigned") <- unassigned
  out
}

# ---- helpers ----------------------------------------------------------------

pick <- function(df, species, iter, source = NULL, years = NULL) {
  keep <- df$species == species & df$iter %in% iter
  if (!is.null(source)) keep <- keep & df$source %in% source
  if (!is.null(years)) keep <- keep & df$year %in% years
  df[keep, , drop = FALSE]
}

lfd_bin_mm <- function(data, species) {
  l <- sort(unique(data$lfd$len_mm[data$lfd$species == species]))
  if (length(l) < 2) return(10)
  min(diff(l))
}
