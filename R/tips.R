# Short explanations of settings, shown as tooltips in the app and usable in
# documentation. One sentence or two each; keep them plain.

floss_tips_list <- list(
  # fishery
  fleets = "Fleets share effort across all species they catch.",
  effort_trend = "Nominal effort over time: start level, peak, later level and the years of the changes.",
  sigmaE = "Year-to-year noise on effort, the same for every species caught by the fleet.",
  creep = "Annual increase in fishing efficiency. It raises F but is not visible in the nominal effort index.",
  # species
  linf = "Asymptotic length (cm) of the von Bertalanffy growth curve.",
  k = "Growth rate (per year): how fast fish approach Linf.",
  M = "Natural mortality (per year).",
  steepness = "Share of unfished recruitment produced at 20% of unfished spawning biomass.",
  sigR = "Recruitment variability (lognormal SD).",
  msy = "Target MSY (t); the stock is scaled so that its equilibrium MSY equals this value.",
  selectivity = "Probability that a fish of a given length is caught by the gear: logistic (flat-topped) or dome-shaped.",
  relq = "Relative catchability: how strongly the fleet fishes this species per unit effort.",
  # operating model
  nsim = "Number of iterations (simulated histories); each can become one dataset.",
  ar1 = "Autocorrelated recruitment deviations: good and bad years come in runs.",
  # treatments
  treatment = "A controlled change from a given year, run against the baseline with the same random numbers.",
  # sampling
  catch_cv = "Observation error of reported catch.",
  catch_bias = "Under-reporting: share of the true catch that is not reported.",
  cpue_rho = "Autocorrelation of CPUE errors. High values mimic non-random sampling, e.g. a fleet fishing the same grounds.",
  cpue_beta = "Proportionality: 1 = CPUE proportional to biomass; below 1 = hyperstable (CPUE stays up as the stock declines).",
  ess_age = "Stage 1 of length sampling: number of independent clusters (trips, hauls) per year.",
  ess_len = "Stage 2 of length sampling: number of fish measured per year.",
  bio_n = "Fish per year with length, age, weight and maturity recorded.",
  age_scheme = "Random: a share of fish is aged. Length-stratified: at most n otoliths per length bin.",
  age_bias = "Systematic ageing error: 0 (default) = unbiased; -5% = fish read 5% too young on average, +5% = 5% too old.",
  q_step = "A sudden change in catchability from a given year (new vessel, gear or electronics).",
  seasons = "Length and biological samples are spread over seasons; fish grow through the year, so the season of sampling changes the lengths seen.",
  season_kappa = "Opportunistic sampling: low values give very uneven shares from year to year, high values stay close to the shares set.",
  survey_sel = "Survey selectivity relative to Linf, so one survey gear fits species of different size.",
  # truth
  hide_truth = "Hide the true stock status and the answer key, e.g. while presenting to trainees."
)

#' Short explanations of FLOSS settings
#'
#' The texts used for the information tooltips of the app.
#' @param id setting name; `NULL` returns all tips.
#' @return a character string (or a named list of all tips).
#' @export
#' @examples
#' floss_tip("cpue_rho")
floss_tip <- function(id = NULL) {
  if (is.null(id)) return(floss_tips_list)
  floss_tips_list[[id]] %||% ""
}
