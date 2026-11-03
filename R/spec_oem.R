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
#' @param q_step artefact: `c(year, mult)`, the survey catchability changes
#'   by `mult` from `year` (e.g. a new vessel or gear), or `NULL`.
#' @return a list of class `floss_survey`.
#' @export
survey_settings <- function(years = NULL, timing = 0, cv = 0.25, rho = 0,
                            sel = list(type = "logistic", L50 = 0.2, L95 = 0.35),
                            lfd = TRUE, ess_age = 50, ess_len = 500,
                            bio_n = 0, q_step = NULL) {
  structure(list(years = years, timing = timing, cv = cv, rho = rho,
                 sel = sel, lfd = lfd, ess_age = ess_age, ess_len = ess_len,
                 bio_n = bio_n, q_step = q_step), class = "floss_survey")
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
#' @param effort_days fishing days per year of a fleet at relative effort 1
#'   (per fleet): observed effort is reported in fishing days.
#' @param cpue_cv,cpue_rho CPUE CV and residual autocorrelation.
#' @param cpue_creep does CPUE carry effort creep (rising catchability)?
#' @param cpue_beta proportionality exponent, \eqn{I = qB^\beta}: 1 =
#'   proportional, < 1 hyperstable (e.g. range contraction), > 1
#'   hyperdepletion.
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
#' @param catch_bias_from year from which `catch_bias` applies, per fleet
#'   (default: all years); an artefact such as under-reporting that starts
#'   with a new quota.
#' @param cpue_q_step artefact: named list by fleet of `c(year, mult)`, the
#'   fleet's catchability changes by `mult` from `year` (e.g. new
#'   electronics), or `NULL`.
#' @param bio_age_bias artefact: proportional ageing bias, e.g. -0.1 = fish
#'   are aged 10% too young (typical of old fish read from whole otoliths).
#' @param survey [survey_settings()] or `NULL` for no survey.
#' @param seasons seasonal sampling of the fleets' length and biological
#'   samples: a [season_design()] for all fleets, or a named list of designs
#'   by fleet. Default: four seasons, balanced. `season_design(1)` gives one
#'   annual sample at mid-year.
#' @param years optional named list by fleet of `c(first, last)` years,
#'   overriding `tail`.
#' @return a list of class `floss_oem`.
#' @export
#' @examples
#' oem_settings(reef_fishery(), tail = c(traps = 5, gillnet = 10, linefish = 15))
oem_settings <- function(fishery, tail = 10, catch_cv = 0.1, catch_bias = 0,
                         effort_cv = 0.1, effort_days = 1000,
                         cpue_cv = 0.2, cpue_rho = 0.3,
                         cpue_creep = TRUE, cpue_beta = 1, lfd_ess_age = 30,
                         lfd_ess_len = 300, lfd_bin = NULL, lfd_cv = 0.1,
                         bio_n = 100, bio_p_aged = 0.5,
                         bio_age_scheme = c("random", "stratified"),
                         bio_n_per_bin = 5, age_cv = 0.1,
                         wt_cv = 0.1, mat_L95 = 1.15, catch_bias_from = NULL,
                         cpue_q_step = NULL, bio_age_bias = 0,
                         survey = survey_settings(), years = NULL,
                         seasons = season_design(4)) {
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
  if (is.numeric(seasons)) seasons <- season_design(seasons)
  if (!inherits(seasons, "floss_seasons")) {
    if (is.null(names(seasons)) || !all(names(seasons) %in% fl))
      stop("'seasons' must be a season_design() or a list of them named by fleet")
    seasons <- lapply(stats::setNames(fl, fl), function(f) seasons[[f]] %||% season_design(4))
  }
  if (!is.null(survey) && is.null(survey$years))
    survey$years <- as.integer(c(last - 24, last))
  structure(list(
    version = floss_spec_version(), fleets = fl,
    years = yrs,
    catch = list(cv = pf(catch_cv, "catch_cv"), bias = pf(catch_bias, "catch_bias"),
                 bias_from = if (!is.null(catch_bias_from))
                   pf(catch_bias_from, "catch_bias_from")),
    effort = list(cv = pf(effort_cv, "effort_cv"), days = pf(effort_days, "effort_days")),
    cpue = list(cv = pf(cpue_cv, "cpue_cv"), rho = pf(cpue_rho, "cpue_rho"),
                creep = cpue_creep, beta = pf(cpue_beta, "cpue_beta"),
                q_step = cpue_q_step),
    lfd = list(ess_age = pf(lfd_ess_age, "lfd_ess_age"),
               ess_len = pf(lfd_ess_len, "lfd_ess_len"),
               bin = lfd_bin, cv = lfd_cv),
    bio = list(n = pf(bio_n, "bio_n"), p_aged = bio_p_aged,
               age_scheme = bio_age_scheme, n_per_bin = bio_n_per_bin,
               age_cv = age_cv, age_bias = bio_age_bias,
               wt_cv = wt_cv, mat_L95 = mat_L95),
    survey = survey, seasons = seasons), class = "floss_oem")
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
  for (f in x$fleets) {
    d <- fleet_seasons(x, f)
    cat("seasons", f, ":", d$n, if (d$n > 1) paste0("(", paste(sprintf("%.2f", d$share), collapse = "/"),
        if (length(d$shift)) "; shifts", if (is.finite(d$kappa)) paste0("; kappa ", d$kappa), ")"), "\n")
  }
  invisible(x)
}

#' Artefacts planted in the sampling settings and the fishery
#'
#' Lists every setting that makes the data depart from the stock in a way a
#' careful analyst should detect (the "warning flags"): effort creep carried
#' by CPUE, non-proportional CPUE, strongly autocorrelated residuals,
#' catch bias, catchability steps, ageing bias and very small length
#' samples, plus treatments applied to the operating model. Written to the
#' answer key by [export_data()].
#'
#' @param oem [oem_settings()].
#' @param fishery the [fishery()].
#' @param treatment optional [treatment()] of the operating model.
#' @return a data frame with `artefact`, `where`, `detail` (empty if none).
#' @export
planted_artefacts <- function(oem, fishery = NULL, treatment = NULL) {
  rows <- list()
  add <- function(a, w, d) rows[[length(rows) + 1]] <<- data.frame(artefact = a, where = w, detail = d)
  for (f in oem$fleets) {
    if (!is.null(fishery) && fishery$creep[[f]] > 0 && isTRUE(oem$cpue$creep))
      add("effort creep in CPUE", f, sprintf("%.1f%% per year from %d", 100 * fishery$creep[[f]],
                                            fishery$creep_start[[f]]))
    if ((oem$cpue$beta[[f]] %||% 1) != 1)
      add("non-proportional CPUE", f, sprintf("beta = %.2f", oem$cpue$beta[[f]]))
    if (oem$cpue$rho[[f]] >= 0.5)
      add("autocorrelated CPUE residuals", f, sprintf("rho = %.2f", oem$cpue$rho[[f]]))
    if (oem$catch$bias[[f]] != 0)
      add("catch bias", f, sprintf("%.0f%% under-reported%s", 100 * oem$catch$bias[[f]],
          if (!is.null(oem$catch$bias_from)) paste(" from", oem$catch$bias_from[[f]]) else ""))
    qs <- oem$cpue$q_step[[f]]
    if (!is.null(qs)) add("CPUE catchability step", f, sprintf("x %.2f from %d", qs[2], qs[1]))
    sd <- fleet_seasons(oem, f)
    if (sd$n > 1 && !is_balanced(sd))
      add("imbalanced seasonal sampling", f, paste0("shares ", paste(sprintf("%.2f", sd$share), collapse = "/"),
          if (length(sd$shift)) paste0("; from ", vapply(sd$shift, function(z) paste0(z$from, ": ",
            paste(sprintf("%.2f", z$share), collapse = "/")), ""), collapse = ""),
          if (is.finite(sd$kappa)) paste0("; random year-to-year, kappa ", sd$kappa)))
    if (!is.null(sd$n_mult) && sd$n_mult$mult != 1)
      add("change in sampling intensity", f, sprintf("sample size x %.2f from %d%s", sd$n_mult$mult,
          sd$n_mult$from, if (sd$n_mult$over > 0) paste0(" over ", sd$n_mult$over, " years") else ""))
    if (sd$p_missing > 0)
      add("missing sampling years", f, sprintf("each year unsampled with probability %.2f", sd$p_missing))
    if (oem$lfd$ess_age[[f]] < 15)
      add("very small length samples", f, sprintf("%d clusters", oem$lfd$ess_age[[f]]))
  }
  if (!is.null(oem$survey$q_step))
    add("survey catchability step", "survey", sprintf("x %.2f from %d",
        oem$survey$q_step[2], oem$survey$q_step[1]))
  if ((oem$bio$age_bias %||% 0) != 0)
    add("ageing bias", "biological samples", sprintf("%+.0f%% of true age", 100 * oem$bio$age_bias))
  if (!is.null(treatment)) add("treatment of the stock", "operating model", treatment$name)
  if (!length(rows)) return(data.frame(artefact = character(), where = character(),
                                       detail = character()))
  do.call(rbind, rows)
}
