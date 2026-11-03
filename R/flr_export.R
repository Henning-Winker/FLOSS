# FLR objects from FLOSS data, ready for JABBA3 (catch, indices, effort) and
# FLicc (length frequencies by gear and life-history parameters).

#' FLR objects for JABBA3 and FLicc
#'
#' Converts the observations of one species and iteration into FLR objects:
#'
#' * `catch`: `FLQuants`, observed catch (t) by fleet (age `"all"`), and
#'   `catch_total` (`FLQuant`, all fleets), as `jb3_data()` of JABBA3 takes
#'   them;
#' * `indices`: `FLIndices` of `FLIndexBiomass`, one per fleet CPUE and the
#'   survey, with the log-scale standard error (from the sampling CV) in
#'   `index.var`, as JABBA3 expects;
#' * `effort`: `FLQuants` of effort by fleet (fishing days), for JABBA3's
#'   `auxiliary` data;
#' * `lfd`: `FLQuants`, one per gear (fleets and survey), length (lower bin
#'   limit, cm) x year, counts summed over seasons, as FLicc's `fiticc()` and
#'   `stocklen()` take them; `lfd_raised` the same raised to the catch
#'   ([lfd_aggregate()]);
#' * `lhpar`: `FLPar` with `linf`, `k`, `t0`, `M`, `L50`, `a`, `b`, `tmax`
#'   and `s` (steepness), with `a` converted so that weight is in kg (FLicc
#'   convention; FLOSS uses g);
#' * `catch_by_gear`: share of the catch by fleet over the years with length
#'   samples (FLicc's `catch_by_gear`).
#'
#' @param data a `floss_data`.
#' @param species species name; default the first. `NULL` with
#'   `all = TRUE` returns a list over species.
#' @param iter one iteration.
#' @param all return a named list with one element per species.
#' @return a list of class `floss_flr` (or a list of them).
#' @export
#' @examples
#' \donttest{
#' x <- simulate_fishery("lutjanus_bohar", nsim = 2, iters = 1)
#' f <- as_flr(x$data)
#' f$indices
#' # JABBA3: jb3_data(f$catch, f$indices, auxiliary = f$effort)
#' # FLicc:  stocklen(f$lfd[c("gillnet", "linefish")], f$lhpar)
#' }
as_flr <- function(data, species = NULL, iter = 1, all = FALSE) {
  if (all) {
    sps <- unique(data$catch$species)
    return(stats::setNames(lapply(sps, function(s) as_flr(data, s, iter)), sps))
  }
  species <- species %||% unique(data$catch$species)[1]
  if (length(iter) != 1) stop("as_flr() works on one iteration")
  oem <- data$settings$oem
  fl <- data$settings$fishery$fleets %||% unique(data$catch$fleet)
  yq <- function(v, yrs, units) FLQuant(as.numeric(v), dimnames = list(age = "all", year = yrs), units = units)
  ct <- data$catch[data$catch$species == species & data$catch$iter == iter, ]
  catch <- FLQuants(stats::setNames(lapply(fl, function(f) {
    d <- ct[ct$fleet == f, ]; d <- d[order(d$year), ]
    yq(d$catch_t, d$year, "t")
  }), fl))
  tot <- stats::aggregate(catch_t ~ year, ct, sum)
  ef <- data$effort[data$effort$iter == iter, ]
  effort <- FLQuants(stats::setNames(lapply(fl, function(f) {
    d <- ef[ef$fleet == f, ]; d <- d[order(d$year), ]
    yq(d$effort, d$year, "days")
  }), fl))
  # indices with log-scale SE from the sampling CV
  idx <- rbind(data$cpue, data$survey)
  idx <- idx[idx$species == species & idx$iter == iter, ]
  cv_of <- function(s) if (s == "survey") oem$survey$cv %||% 0.25 else oem$cpue$cv[[s]] %||% 0.2
  srcs <- intersect(c(fl, "survey"), unique(idx$source))
  indices <- FLIndices(stats::setNames(lapply(srcs, function(s) {
    d <- idx[idx$source == s, ]; d <- d[order(d$year), ]
    fi <- FLIndexBiomass(index = yq(d$index, d$year, "NA"))
    index.var(fi) <- yq(rep(sqrt(log(1 + cv_of(s)^2)), nrow(d)), d$year, "NA")
    name(fi) <- s
    fi
  }), srcs))
  # length frequencies by gear on one length grid
  bw <- lfd_bin_mm(data, species)
  lf_flqs <- function(l) {
    if (is.null(l) || !nrow(l)) return(NULL)
    grid <- seq(min(l$len_mm), max(l$len_mm), by = bw)
    gs <- intersect(c(fl, "survey"), unique(l$source))
    FLQuants(stats::setNames(lapply(gs, function(g) {
      d <- l[l$source == g, ]
      yrs <- sort(unique(d$year))
      m <- matrix(0, length(grid), length(yrs), dimnames = list(grid / 10, yrs))
      m[cbind(match(d$len_mm, grid), match(d$year, yrs))] <- d$n
      FLQuant(m, dimnames = list(len = as.character(grid / 10), year = yrs), units = "cm")
    }), gs))
  }
  l_s <- pick(lfd_annual(data$lfd), species, iter)
  l_r <- pick(lfd_aggregate(data, species, "fleet", "catch"), species, iter)
  lh <- data$lhpars[data$lhpars$species == species, ]
  lhpar <- FLPar(linf = lh$linf_cm, k = lh$k, t0 = lh$t0, M = lh$M, L50 = lh$lm_cm,
                 a = lh$a / 1000, b = lh$b, tmax = lh$tmax, s = lh$steepness)
  ly <- unique(l_s$year[l_s$source != "survey"])
  sh <- stats::aggregate(catch_t ~ fleet, ct[ct$year %in% ly, ], sum)
  share <- stats::setNames(sh$catch_t / sum(sh$catch_t), sh$fleet)[fl]
  structure(list(species = species, iter = iter, catch = catch,
                 catch_total = yq(tot$catch_t, tot$year, "t"),
                 indices = indices, effort = effort,
                 lfd = lf_flqs(l_s), lfd_raised = lf_flqs(l_r),
                 lhpar = lhpar, catch_by_gear = share[!is.na(share)]),
            class = "floss_flr")
}

#' @exportS3Method base::print
print.floss_flr <- function(x, ...) {
  cat("FLOSS data as FLR objects:", x$species, "- iteration", x$iter, "\n")
  cat("  catch (FLQuants, t):", paste(names(x$catch), collapse = ", "), "\n")
  cat("  indices (FLIndices):", paste(names(x$indices), collapse = ", "), "\n")
  cat("  effort (FLQuants, days):", paste(names(x$effort), collapse = ", "), "\n")
  cat("  lfd (FLQuants, len x year):", paste(names(x$lfd), collapse = ", "), "\n")
  cat("  lhpar:", paste(dimnames(x$lhpar)$params, collapse = ", "), "\n")
  cat("  JABBA3: jb3_data(x$catch, x$indices, auxiliary = x$effort)\n")
  cat("  FLicc:  stocklen(x$lfd, x$lhpar); catch_by_gear = x$catch_by_gear\n")
  invisible(x)
}
