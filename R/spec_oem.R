# Sampling (observation error model) settings.

#' Survey settings
#'
#' A fishery-independent survey that samples every species. Because one survey
#' gear meets species of very different size, its selectivity is given
#' **relative to each species' L-infinity** (e.g. `L50 = 0.2` means 50%
#' selected at 20% of Linf).
#'
#' @param years first and last survey year (default: set by
#'   [oem_settings()] to the last 25 years).
#' @param timing fraction of the year when the survey takes place (0 = 1 Jan).
#' @param cv,rho observation CV and lag-1 autocorrelation of the index
#'   residuals.
#' @param sel selectivity relative to Linf, as for [selectivity()].
#' @param lfd,ess_age,ess_len take two-stage length samples from the survey
#'   (number of stations and of fish measured).
#' @param bio_n fish per year for biological samples from the survey (0 =
#'   none).
#' @return a list of class `floss_survey`.
#' @export
survey_settings <- function(years = NULL, timing = 0, cv = 0.25, rho = 0,
                            sel = list(type = "logistic", L50 = 0.2, L95 = 0.35),
                            lfd = TRUE, ess_age = 50, ess_len = 500,
                            bio_n = 0) {
  structure(list(years = years, timing = timing, cv = cv, rho = rho,
                 sel = sel, lfd = lfd, ess_age = ess_age, ess_len = ess_len,
                 bio_n = bio_n), class = "floss_survey")
}

#' Sampling settings
#'
#' Collects every setting of the observation model. Arguments that can differ
#' by fleet take one value or a vector named by fleet. Catch and effort cover
#' all years of the operating model; CPUE, length frequencies and biological
#' samples cover the last `tail` years of each fleet (or `years`).
#'
#' * **Catch** by fleet: true catch x (1 - bias) x lognormal error.
#' * **Effort** index by fleet: nominal effort x lognormal error (effort creep
#'   is not visible).
#' * **CPUE** by fleet: vulnerable biomass x efficiency (if `cpue_creep`) with
#'   lognormal AR1 residuals,
#'   \eqn{\epsilon_y = \rho\,\epsilon_{y-1} + \sqrt{1-\rho^2}\,\delta_y}. Use
#'   \eqn{\rho > 0} for non-random (clustered) sampling of the fleet.
#' * **Length frequencies**: the two-stage FLRef protocol
#'   ([FLRef::lfd.sim()]): `lfd_ess_age` clusters drawn from the fleet's
#'   catch-at-age, then `lfd_ess_len` fish measured through the gear's
#'   length key.
#' * **Biological samples**: `bio_n` fish per fleet and year with length
#'   (mm), age (with ageing error; aged at random or length-stratified),
#'   weight (g, with error) and maturity (0/1).
#'
#' @param fishery a [fishery()] (or a `floss_om` / `floss_oms`).
#' @param tail number of most recent years with CPUE, length and biological
#'   samples, per fleet.
#' @param catch_cv,catch_bias catch CV and under-reporting (0.1 = 10% under).
#' @param effort_cv effort-index CV.
#' @param cpue_cv,cpue_rho CPUE CV and residual autocorrelation.
#' @param cpue_creep does CPUE carry effort creep (hyperstability)?
#' @param lfd_ess_age,lfd_ess_len stage-1 (clusters) and stage-2 (fish
#'   measured) sample sizes.
#' @param lfd_bin length bin (cm); `NULL` picks about Linf/40 per species.
#' @param lfd_cv CV of length-at-age.
#' @param bio_n,bio_p_aged fish per year and share aged (random ageing).
#' @param bio_age_scheme `"random"` (a share `bio_p_aged` of fish aged) or
#'   `"stratified"` (length-stratified subsample: at most `bio_n_per_bin`
#'   fish aged per length bin, source and year, as [FLRef::alk.sample()]).
#' @param bio_n_per_bin otoliths per length bin for the stratified scheme.
#' @param age_cv,wt_cv ageing-error CV and weight CV.
#' @param mat_L95 length at 95% maturity relative to `lm` (L50).
#' @param survey [survey_settings()] or `NULL` for no survey.
#' @param years optional named list by fleet of `c(first, last)` years,
#'   overriding `tail`.
#' @return a list of class `floss_oem`.
#' @export
#' @examples
#' oem_settings(reef_fishery(), tail = c(traps = 5, gillnet = 10, linefish = 15))
oem_settings <- function(fishery, tail = 10, catch_cv = 0.1, catch_bias = 0,
                         effort_cv = 0.1, cpue_cv = 0.2, cpue_rho = 0.3,
                         cpue_creep = TRUE, lfd_ess_age = 30,
                         lfd_ess_len = 300, lfd_bin = NULL, lfd_cv = 0.1,
                         bio_n = 100, bio_p_aged = 0.5,
                         bio_age_scheme = c("random", "stratified"),
                         bio_n_per_bin = 5, age_cv = 0.1,
                         wt_cv = 0.1, mat_L95 = 1.15,
                         survey = survey_settings(), years = NULL) {
  bio_age_scheme <- match.arg(bio_age_scheme)
  if (inherits(fishery, "floss_om")) fishery <- fishery$fishery
  if (inherits(fishery, "floss_oms")) fishery <- fishery[[1]]$fishery
  fl <- fishery$fleets
  last <- fishery$years[2]
  pf <- function(x, what) {
    if (length(x) == 1 && is.null(names(x))) x <- stats::setNames(rep(x, length(fl)), fl)
    if (is.null(names(x))) names(x) <- fl
    miss <- setdiff(fl, names(x))
    if (length(miss)) stop("'", what, "' has no value for ", toString(miss))
    as.list(x[fl])
  }
  tl <- pf(tail, "tail")
  yrs <- lapply(fl, function(f) {
    if (!is.null(years[[f]])) return(as.integer(years[[f]]))
    as.integer(c(last - tl[[f]] + 1, last))
  })
  names(yrs) <- fl
  if (!is.null(survey) && is.null(survey$years))
    survey$years <- as.integer(c(last - 24, last))
  structure(list(
    version = floss_spec_version(), fleets = fl,
    years = yrs,
    catch = list(cv = pf(catch_cv, "catch_cv"), bias = pf(catch_bias, "catch_bias")),
    effort = list(cv = pf(effort_cv, "effort_cv")),
    cpue = list(cv = pf(cpue_cv, "cpue_cv"), rho = pf(cpue_rho, "cpue_rho"),
                creep = cpue_creep),
    lfd = list(ess_age = pf(lfd_ess_age, "lfd_ess_age"),
               ess_len = pf(lfd_ess_len, "lfd_ess_len"),
               bin = lfd_bin, cv = lfd_cv),
    bio = list(n = pf(bio_n, "bio_n"), p_aged = bio_p_aged,
               age_scheme = bio_age_scheme, n_per_bin = bio_n_per_bin,
               age_cv = age_cv,
               wt_cv = wt_cv, mat_L95 = mat_L95),
    survey = survey), class = "floss_oem")
}

#' @exportS3Method base::print
print.floss_oem <- function(x, ...) {
  cat("FLOSS sampling settings\n")
  tab <- data.frame(
    fleet = x$fleets,
    years = vapply(x$years, paste, "", collapse = "-"),
    catch_cv = unlist(x$catch$cv), bias = unlist(x$catch$bias),
    effort_cv = unlist(x$effort$cv), cpue_cv = unlist(x$cpue$cv),
    cpue_rho = unlist(x$cpue$rho), ess_age = unlist(x$lfd$ess_age),
    ess_len = unlist(x$lfd$ess_len), bio_n = unlist(x$bio$n), row.names = NULL)
  print(tab, row.names = FALSE)
  if (!is.null(x$survey))
    cat("survey:", paste(x$survey$years, collapse = "-"), " timing",
        x$survey$timing, " cv", x$survey$cv, " rho", x$survey$rho,
        if (isTRUE(x$survey$lfd)) " + length samples", "\n")
  invisible(x)
}
