# =====================================================================================================================
# FLOSS user guide as an R script
# ---------------------------------------------------------------------------------------------------------------------
# Four reef species and three fleets, from the fishery to saved datasets: fishery and effort,
# species, operating models, sampling (seasons, length and biological samples), observations,
# naive aggregation for teaching, true status, export (csv, Excel, FLR objects for JABBA3 and
# FLicc), treatments and an effort-creep teaching scenario.
#
# Generated from inst/guide/user-guide.Rmd with knitr::purl(documentation = 2) (dev/make_guide_scripts.R): the guide text
# is kept as #' comments, section headings as "#' #" lines, and every code chunk is marked
# "## ---- <chunk name> ----". Chunks the guide shows but does not run (eval = FALSE) are commented out with "## ".
#
# Requirements: FLOSS and its FLR dependencies; writexl for Excel export; FLicc and JABBA3 only for the commented-out examples.
# Install: install.packages(c("FLCore", "ggplotFL", "FLBRP", "FLasher", "FLife", "FLRef"),
#                           repos = "https://flr.r-universe.dev")
#          remotes::install_github("Henning-Winker/FLOSS")
# Run section by section (RStudio: Ctrl+Enter, or Ctrl+Alt+T for the current chunk). Objects are reused in later
# sections, so run from the top. The rendered guide: https://henning-winker.github.io/FLOSS/
# =====================================================================================================================

#' 
#' This guide runs the whole FLOSS workflow for the multi-gear reef example:
#' four species (*Siganus sutor*, *Lethrinus mahsena*, *Lutjanus bohar*,
#' *Carcharhinus melanopterus*) caught by three fleets (traps, gillnet,
#' linefish). It ends with datasets saved to disk, one per iteration, with a
#' separate answer key: the material a trainer hands out in a stock-assessment
#' course.
#' 
#' ```
#' 1 Fishery -> 2 Species -> 3 Operating models -> 4 Sampling -> 5 Data -> 6 Observations -> 7 Truth -> 8 Save
#' ```
#' 
#' NOTE:
#' **Where to read more.** `vignette("floss-00-how-it-works")` describes the
#' model, the sampling and what to expect in the data;
#' `vignette("floss-03-teaching")` has the teaching cards (predict the
#' response, spot the warning flags). This guide is the step-by-step workflow.
#' 
## ----intro--------------------------------------------------------------------
library(FLOSS)
# keep species saved in this guide out of your own library
options(FLOSS.user_dir = file.path(tempdir(), "floss_user"))

#' 
#' # 1 Fishery
#' 
#' A fishery is a set of fleets with nominal effort trends over a common time
#' horizon. It is species-independent: every species sees the same effort, so
#' species caught by the same fleet share effort anomalies.
#' 
## ----1-fishery----------------------------------------------------------------
fy <- load_fishery("reef_3fleets")
fy

#' 
## ----1-fishery-2, fig.height = 6----------------------------------------------
plot_effort(fy, nsim = 50)

#' 
#' NOTE:
#' **Shared effort.** Fishing mortality of species *s* by fleet *f* is
#' `relq[s, f] × Fmsy[s] × efficiency[f, y] × effort[f, y]`. Effort noise is
#' common to all species of a fleet; effort creep raises efficiency but is not
#' visible in the effort index.
#' 
#' Your own fishery is defined with `fishery()` and `effort_trend()`; for
#' example, two fleets with 1% annual effort creep in the second:
#' 
## ----1-fishery-3, eval = FALSE------------------------------------------------
## my_fy <- fishery(
##   fleets = c("handline", "net"), years = 1980:2025,
##   trends = list(handline = effort_trend(fhi = 1.5, flo = 1, breaks = c(0.5, 0.8)),
##                 net      = effort_trend(fhi = 0.8, flo = 1.2, breaks = c(0.3, 0.6))),
##   sigmaE = 0.15, creep = c(handline = 0, net = 0.01))
## write_spec(my_fy, "my_fishery.yaml")

#' 
#' # 2 Species
#' 
#' FLOSS ships the four reef species as read-only examples. Each holds
#' life-history parameters (FishLife), a target MSY, selectivity by fleet and
#' relative catchability by fleet.
#' 
## ----2-species----------------------------------------------------------------
species_library()[, c("id", "common_name", "msy", "source")]
sp <- load_species("lutjanus_bohar")
sp

#' 
## ----2-species-2, fig.height = 3.5--------------------------------------------
FLRef::plot_sel_la(species_sel(sp), age_by = 2, len_by = 10)

#' 
#' A new species is easiest made from a copy. Here a hypothetical humpback
#' snapper borrows the selectivity of *L. bohar*:
#' 
## ----2-species-3--------------------------------------------------------------
gibbus <- species(
  "Lutjanus gibbus", common_name = "Humpback red snapper",
  lh = modifyList(sp$lh, list(linf = 45, k = 0.3, lm = 25, tm = 2.5, tmax = 18,
                              M = 0.35)),
  msy = 400, sel = sp$sel, relq = c(traps = 0.4, gillnet = 0.6, linefish = 0.8))
save_species(gibbus)
species_library()[, c("id", "msy", "source")]

#' 
#' The rest of the guide uses the four example species.
#' 
#' # 3 Operating models
#' 
#' `build_oms()` calibrates each species to its MSY, applies fleet F (effort ×
#' relative catchability × Fmsy × selectivity) and projects the stock with
#' stochastic recruitment. Use 100 or more iterations in practice; 20 keeps this
#' guide quick.
#' 
## ----3-operating-models-------------------------------------------------------
ids <- c("siganus_sutor", "lethrinus_mahsena", "lutjanus_bohar",
         "carcharhinus_melanopterus")
oms <- build_oms(ids, fy, nsim = 20)
oms$lutjanus_bohar

#' 
## ----3-operating-models-2, fig.height = 6-------------------------------------
plot_om(oms)

#' 
## ----3-operating-models-3, fig.height = 3-------------------------------------
plot_fleets(oms$lutjanus_bohar)

#' 
#' # 4 Sampling settings
#' 
#' `oem_settings()` holds every observation setting; `reef_oem()` gives the
#' settings of this example: CPUE, length frequencies and biological samples for
#' the most recent **5 years (traps), 10 years (gillnet) and 15 years
#' (linefish)**, catch and effort for all years, and a January survey for the
#' last 25 years.
#' 
## ----4-sampling-settings------------------------------------------------------
oem <- reef_oem()
oem

#' 
#' NOTE:
#' **Why sampling settings matter.** The same stock looks different through
#' different sampling: few clusters make length samples erratic, autocorrelated
#' CPUE errors make trends look more certain than they are, and the years
#' sampled decide whether a change in the stock is visible at all.
#' 
#' Three choices matter most:
#' 
#' * **CPUE residual autocorrelation** (`cpue_rho`). A fleet that keeps fishing
#'   the same grounds does not sample the stock at random, so its CPUE errors
#'   persist from year to year:
#'   \eqn{\epsilon_y = \rho\,\epsilon_{y-1} + \sqrt{1-\rho^2}\,\delta_y}. The
#'   linefish fleet uses \eqn{\rho = 0.6}.
#' * **Two-stage length sampling** (`lfd_ess_age`, `lfd_ess_len`), the FLRef
#'   protocol of `FLRef::lfd.sim()`: clusters (trips, hauls) are drawn from the
#'   fleet's catch-at-age, then fish are measured through the gear's length
#'   key. The effective sample size lies between the two.
#' * **Ageing scheme** of the biological samples (`bio_age_scheme`): a random
#'   share of fish (`"random"`, `bio_p_aged`) or a length-stratified subsample
#'   of at most `bio_n_per_bin` otoliths per length bin (`"stratified"`).
#' 
#' Any setting can be changed, by fleet where it makes sense:
#' 
## ----4-sampling-settings-2----------------------------------------------------
oem_strat <- reef_oem(bio_age_scheme = "stratified", bio_n_per_bin = 5,
                      catch_bias = c(traps = 0, gillnet = 0, linefish = 0.2))

#' 
#' # 5 Generate the data
#' 
#' `sample_data()` samples every data stream for the chosen iterations. Each
#' iteration becomes one dataset (e.g. one per group of trainees).
#' 
## ----5-generate-the-data------------------------------------------------------
dat <- sample_data(oms, oem, iters = 1:3, seed = 2026)
dat

#' 
#' The result is a list of plain data frames:
#' 
## ----5-generate-the-data-2----------------------------------------------------
head(dat$catch, 3)
head(dat$cpue, 3)
head(dat$lfd, 3)
head(dat$bio, 3)

#' 
#' # 6 Look at the observations
#' 
#' All plots show observations only; nothing is fitted.
#' 
#' ## Which data exist for which years
#' 
## ----which-data-exist-for, fig.height = 5-------------------------------------
plot_availability(dat, "Lutjanus bohar")

#' 
#' ## Catch, CPUE and survey
#' 
## ----catch-cpue-and-survey----------------------------------------------------
plot_obs(dat, "Lutjanus bohar", what = "catch")

#' 
#' With `truth = TRUE` (trainer view) the true index is added as a dashed line:
#' 
## ----catch-cpue-and-survey-2--------------------------------------------------
plot_obs(dat, "Lutjanus bohar", what = "index", truth = TRUE)

#' 
#' Several datasets can be compared by passing more than one iteration:
#' 
## ----catch-cpue-and-survey-3--------------------------------------------------
plot_obs(dat, "Lutjanus bohar", iter = 1:3)

#' 
#' ## Length frequencies
#' 
#' The classic view: one length-frequency histogram per year and source. The
#' dotted lines mark the mean length of each age class in mid-year, so modes
#' can be matched to year classes.
#' 
## ----length-frequencies, fig.height = 7, fig.width = 8------------------------
plot_lfd(dat, "Lutjanus bohar", type = "bars", source = c("gillnet", "linefish", "survey"),
         years = 2018:2025, cohorts = TRUE)

#' 
#' The ridge view shows each year's length distribution along the time axis
#' with von Bertalanffy cohort lines: real modes move along them.
#' 
## ----length-frequencies-2, fig.height = 7-------------------------------------
plot_lfd(dat, "Lutjanus bohar", type = "ridge", source = c("gillnet", "linefish"), cohorts = TRUE)

#' 
#' The bubble view shows all years at once. A strong year class appears as a
#' ridge of large bubbles moving to the right along the growth lines.
#' 
## ----length-frequencies-2b, fig.height = 5.5----------------------------------
plot_lfd(dat, "Lutjanus bohar", type = "bubble", cohorts = TRUE)

#' 
#' NOTE:
#' **Seasonal sampling.** The samples of each year are spread over four
#' seasons (the default, `season_design()`). Fish grow through the year, so
#' later seasons have larger fish. `plot_seasons()` shows the sample sizes and
#' the mean length by season; `plot_lfd(season = 1)` shows one season only.
#' 
## ----length-frequencies-3, fig.height = 5.5-----------------------------------
plot_seasons(dat, "Lutjanus bohar")

#' 
#' With `by_season = TRUE` the length samples of each season are shown
#' separately; here an imbalanced design where sampling shifts to the second
#' half of the year from 2020 (set in `oem_settings(seasons = season_design(...))`,
#' or in the app under Sampling, Seasons):
#' 
## ----seasonal-lfd, fig.height = 6.5, fig.width = 8----------------------------
oem_shift <- reef_oem(seasons = season_design(share = c(0.6, 0.3, 0.05, 0.05),
                                              shift = list(from = 2020, share = c(0.05, 0.05, 0.3, 0.6))))
dat_shift <- sample_data(oms[["lutjanus_bohar"]], oem_shift, iters = 1)
plot_lfd(dat_shift, "Lutjanus bohar", type = "ridge", source = "linefish",
         years = 2016:2025, cohorts = TRUE, by_season = TRUE)
plot_lfd(dat_shift, "Lutjanus bohar", type = "bars", source = "linefish",
         years = 2018:2021, by_season = TRUE)

#' 
#' ## Biological samples
#' 
## ----biological-samples, fig.height = 4---------------------------------------
plot_bio(dat, "Lutjanus bohar", what = "age_length")
plot_bio(dat, "Lutjanus bohar", what = "length_weight")
plot_bio(dat, "Lutjanus bohar", what = "maturity")

#' 
#' ## Age-length key and catch-at-age
#' 
#' The observed age-length key counts aged fish by length bin and age:
#' 
## ----age-length-key-and-------------------------------------------------------
alk <- alk_obs(dat, "Lutjanus bohar", source = "linefish")
alk[1:8, 1:8]

#' 
## ----age-length-key-and-2, fig.height = 4-------------------------------------
plot_bio(dat, "Lutjanus bohar", what = "alk", source = "linefish")

#' 
#' Applying the key to the length samples gives catch-at-age. Length bins
#' without aged fish cannot be assigned and are reported as unassigned. With
#' `raise = "catch"` the numbers are raised to the observed catch weight.
#' 
## ----age-length-key-and-3-----------------------------------------------------
caa <- caa_obs(dat, "Lutjanus bohar", source = "linefish", raise = "catch")
round(attr(caa, "unassigned"), 3)

#' 
## ----age-length-key-and-4, fig.height = 4-------------------------------------
plot_caa(caa)

#' 
#' WARNING FLAG:
#' **Before trusting a single indicator**, check that cohorts move through the
#' length data at the growth rate and that the indices broadly agree.
#' `plot_indicators()` puts all series on one page.
#' 
## ----age-length-key-and-5, fig.height = 6-------------------------------------
plot_indicators(dat, "Lutjanus bohar")

#' 
#' # 7 True status (trainers)
#' 
#' The answer key summarised: B/Bmsy, F/Fmsy and catch/MSY in the last year and
#' the Kobe status, for each species and dataset, or as medians over datasets.
#' 
## ----7-true-status-trainers---------------------------------------------------
status_table(dat)
status_table(dat, summary = TRUE)

#' 
## ----7-true-status-trainers-2, fig.height = 6.5, fig.width = 9----------------
plot_truth(dat)

#' 
## ----7-true-status-trainers-3, fig.height = 4.5, fig.width = 6.5--------------
plot_kobe(dat)

#' 
#' The true time series, reference points and status table can be saved on
#' their own, e.g. for a trainer's slides:
#' 
## ----7-true-status-trainers-4-------------------------------------------------
tdir <- export_truth(dat, file.path(tempdir(), "reef_truth"), format = "csv")
list.files(tdir)

#' 
#' ## Naive aggregation (teaching)
#' 
#' An analyst often has catch and effort for the whole fishery, not by fleet,
#' and sometimes for all species together. FLOSS reports effort in fishing
#' days, so catch, effort and nominal CPUE (kg per fishing day) can be shown
#' by fleet, summed over fleets, or summed over fleets and species.
#' 
## ----naive-aggregation-teaching, fig.height = 4-------------------------------
plot_aggregate(dat, "Lutjanus bohar", what = "cpue", by = "fleet", truth = TRUE)

#' 
## ----naive-aggregation-teaching-2, fig.height = 3.5, fig.width = 6------------
plot_aggregate(dat, "Lutjanus bohar", what = "cpue", by = "aggregate", truth = TRUE)
plot_aggregate(dat, what = "cpue", by = "combined")

#' 
#' WARNING FLAG:
#' **Warning flag.** Nominal CPUE summed over fleets changes whenever the fleet
#' mix changes (here linefish effort rises and falls), even if the stock does
#' not; summed over species, the decline of one species can hide behind the
#' others. Compare with the fleet-level series and, as a trainer, the true
#' spawning biomass (grey).
#' 
#' The tables behind these plots are `catch_table()`, `effort_table()` and
#' `nominal_cpue()`.
#' 
#' Length samples can be pooled over fleets as measured, or raised to each
#' fleet's catch first (numbers = sample numbers x catch / sample weight):
#' 
## ----naive-aggregation-teaching-3, fig.height = 5.5, fig.width = 7------------
pooled <- lfd_aggregate(dat, "Lutjanus bohar", by = "aggregate", raise = "sample")
raised <- lfd_aggregate(dat, "Lutjanus bohar", by = "aggregate", raise = "catch")
mean_len <- function(l) tapply(l$len_mm * l$n, l$year, sum) / tapply(l$n, l$year, sum) / 10
round(rbind(pooled = mean_len(pooled[pooled$iter == 1, ]),
            raised = mean_len(raised[raised$iter == 1, ])), 1)
plot_lfd(dat, "Lutjanus bohar", type = "bars", lfd = raised, years = 2020:2025)

#' 
#' Pooled as measured, each fleet counts by its sample size; raised, by its
#' catch. Biological samples can be shown with a panel per fleet, all fleets in
#' one panel, or a panel per season (points coloured by fleet throughout).
#' Length at age rises through the year; as a trainer, `truth = TRUE` adds the
#' true growth curve at the middle of each season:
#' 
## ----naive-aggregation-teaching-4, fig.height = 4-----------------------------
plot_bio(dat, "Lutjanus bohar", what = "age_length", by = "season", truth = TRUE)

#' 
#' # 8 Save the datasets
#' 
#' `export_data()` writes one folder per iteration and the answer key in its own
#' folder; `give = "data"` leaves out the life-history parameters,
#' `format = c("csv", "xlsx")` adds one Excel workbook per iteration, and
#' `zip = TRUE` writes a zip file to hand out.
#' 
## ----8-save-the-datasets------------------------------------------------------
out <- export_data(dat, file.path(tempdir(), "reef_course"), give = "data_lh",
                   format = "csv", answer_key = TRUE)
list.files(out, recursive = TRUE)

#' 
## ----8-save-the-datasets-2, comment = ""--------------------------------------
cat(readLines(file.path(out, "README.txt")), sep = "\n")

#' 
#' Trainees (or you) can read a dataset back and use the same plots:
#' 
## ----8-save-the-datasets-3, fig.height = 4------------------------------------
d1 <- read_dataset(file.path(out, "iter_01"))
names(d1)
plot_bio(d1, "Siganus sutor", what = "age_length")

#' 
#' `settings.yaml` holds the fishery, species and sampling settings and the seed,
#' so the datasets can be re-created exactly. `aggregate = TRUE` adds the
#' aggregated tables (catch and effort summed over fleets and species, nominal
#' CPUE, pooled and raised length samples).
#' 
#' ## FLR objects for JABBA3 and FLicc
#' 
#' `format = "rds"` writes, per iteration, `floss_flr.rds`: one `as_flr()`
#' object per species with FLQuants of catch, effort and length frequencies,
#' FLIndices (with log-scale SEs in `index.var`) and the life-history FLPar.
#' 
## ----flr-objects-for-jabba3---------------------------------------------------
f <- as_flr(dat, "Lutjanus bohar", iter = 1)
f

#' 
## ----flr-objects-for-jabba3-2, eval = FALSE-----------------------------------
## # production model (JABBA3): catch by fleet, all indices, effort as auxiliary data
## jb <- JABBA3::jb3_data(f$catch, f$indices, stock = "bohar", auxiliary = f$effort)
## # length-based (FLicc): length frequencies by gear and life history
## stk <- FLicc::stocklen(f$lfd[c("gillnet", "linefish")], f$lhpar)
## fit <- FLicc::fiticc(f$lfd[c("gillnet", "linefish")], stk,
##                      catch_by_gear = f$catch_by_gear[c("gillnet", "linefish")] /
##                        sum(f$catch_by_gear[c("gillnet", "linefish")]))

#' 
#' # 9 Treatments: one change, same random numbers
#' 
#' A treatment changes one thing in the stock from a given year (effort,
#' recruitment, selectivity, creep, M). `run_treatment()` runs it against the
#' baseline with the same random numbers; see `vignette("floss-03-teaching")`
#' for the teaching cards.
#' 
## ----9-treatments-one-change, fig.height = 5.5--------------------------------
tr <- run_treatment("lutjanus_bohar", fy,
                    treatment("halve linefish effort from 2010",
                              effort = list(fleet = "linefish", from = 2010, mult = 0.5)),
                    nsim = 10)
plot_treatment(tr)

#' 
#' # 10 Everything in one call
#' 
#' `simulate_fishery()` builds the operating models and samples the data in one
#' step:
#' 
## ----10-everything-in-one, eval = FALSE---------------------------------------
## sim <- simulate_fishery(ids, fishery = load_fishery(), oem = reef_oem(),
##                         nsim = 100, iters = 1:10)
## export_data(sim$data, "reef_course", format = c("csv", "xlsx"), zip = TRUE)

#' 
#' # 11 A teaching scenario: effort creep
#' 
#' With 3% annual efficiency increase in the linefish fleet, nominal effort is
#' unchanged but fishing mortality rises. Linefish CPUE, which carries the
#' creep, stays high while the stock declines (hyperstability):
#' 
## ----11-a-teaching-scenario---------------------------------------------------
fy_creep <- reef_fishery(creep = c(traps = 0, gillnet = 0, linefish = 0.03))
om_creep <- build_om(load_species("lutjanus_bohar"), fy_creep, nsim = 20)
dat_creep <- sample_data(om_creep, reef_oem(fishery = fy_creep), iters = 1)
plot_obs(dat_creep, what = "index", truth = TRUE)

#' 
#' Compare `oms$lutjanus_bohar` and `om_creep` with `plot_om()`: the stock is
#' more depleted, but the linefish CPUE barely shows it.
