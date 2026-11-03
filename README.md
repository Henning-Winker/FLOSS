# FLOSS <img src="https://img.shields.io/badge/status-v0.4.1%20prototype-orange" align="right"/>

**FL Operating model & Sampling Simulator**

FLOSS builds multi-fleet, multi-species age-structured operating models and
the data such fisheries produce: catch and effort by fleet, CPUE,
survey indices, length frequencies and biological samples (length, age,
weight, maturity) with realistic sampling error. It is made for **training**
(trainees analyse data from a fishery whose truth is known) and for
**trainers** (generate example datasets in a few clicks, with an answer key).

It builds on [FLR](https://flr-project.org) (`FLCore`, `FLBRP`, `FLasher`,
`FLife`) and [`FLRef`](https://github.com/flr/FLRef).

## Concept

```
1 FISHERY        ->  2 SPECIES           ->  3 OM                 ->  4 SAMPLING      ->  5 EXPORT
fleets, years,       life history, MSY,      calibrate to MSY,        catch, effort,      pick species
nominal effort,      selectivity and         nsim iterations,         CPUE, survey,       and iterations
effort creep         catchability by fleet   truth + refpts           LFD, bio samples    -> csv/xlsx/rds
```

* **Shared effort.** Fleets have nominal effort trends with process noise that
  is common to all species, so species caught by the same fleet share effort
  anomalies. Apical F of species *s* and fleet *f* is
  `F[s,f,y] = relq[s,f] * Fmsy[s] * efficiency[f,y] * effort[f,y]`.
* **Effort creep.** An annual efficiency increase raises F but is not visible
  in the nominal effort index, which produces hyperstable CPUE.
* **Species library.** Four example reef species (*Siganus sutor*,
  *Lethrinus mahsena*, *Lutjanus bohar*, *Carcharhinus melanopterus*); users
  can copy, edit and save their own as plain YAML.
* **Engine first.** Everything runs from R scripts; the Shiny app
  (`run_app()`) only collects inputs and calls the same functions.

## The app

```r
install.packages(c("shiny", "shinydashboard", "commonmark"))
FLOSS::run_app()
```

The app has the layout of the FLicc app: numbered steps in the left panel,
sub-tabs on each page, a *Before you run* checklist with the **Build
operating model** and **Sample data** buttons, and a status pill at the top.

| Step | Sub-tabs |
|---|---|
| Guide | Start here, How it works, Sampling and size data, Warning flags, Truth and export, Settings explained |
| 1 Fishery | Fleets and effort, Effort creep, Effort realisation |
| 2 Species | Library (choose, add, save), Life history, Selectivity, Exploitation |
| 3 Operating model | Settings, Stock trends, Fleets, Equilibrium, Spawners vs vulnerable, Stock and recruitment |
| 4 Treatments | Build, Predict, Reveal (prediction vs what happened, true response, length data) |
| 5 Sampling | Overview, Catch, Effort, CPUE, Survey, Length frequencies, Biological samples |
| 6 Export | Datasets (zip of csv/xlsx with optional answer key), Settings (YAML) |
| 7 True status | Status, Trajectories, Kobe plot, Downloads (only when **Hide truth** is unticked) |

## Installation

```r
install.packages(c("FLCore", "ggplotFL", "FLBRP", "FLasher", "FLife", "FLRef"),
                 repos = c("https://flr.r-universe.dev", "https://cloud.r-project.org"))
remotes::install_github("Henning-Winker/FLOSS")
```

## Quick start

```r
library(FLOSS)

# 1-3: fishery, species and operating models (4 reef species, 3 fleets)
fy  <- load_fishery("reef_3fleets")
oms <- build_oms(species_library()$id, fy, nsim = 100)
plot_om(oms)

# 4-5: sampling settings and data (CPUE/length/bio: last 5, 10, 15 years)
oem <- reef_oem()
dat <- sample_data(oms, oem, iters = 1:10)
plot_availability(dat, "Lutjanus bohar")
plot_obs(dat, "Lutjanus bohar", what = "index", truth = TRUE)
plot_lfd(dat, "Lutjanus bohar")
plot_bio(dat, "Lutjanus bohar", what = "alk")
plot_caa(caa_obs(dat, "Lutjanus bohar", "linefish"))

# 6: save one dataset per iteration + answer key
export_data(dat, "reef_course", format = c("csv", "xlsx"), zip = TRUE)
```

Vignettes:

* `vignette("floss-00-how-it-works")`: the model, sampling, what to expect in size data, warning flags
* `vignette("floss-01-reef-om")`: the reef operating model
* `vignette("floss-02-user-guide")`: step-by-step workflow, four species, three fleets, saved datasets
* `vignette("floss-03-teaching")`: teaching cards (predict the response, spot the warning flags)

### Teaching: treatments with common random numbers

```r
sp <- one_gear_species(load_species("lutjanus_bohar"))   # handline; effort = F/Fmsy
rt <- run_treatment(sp, one_gear_fishery(1),
        list(pulse   = treatment("strong 2012 year class",
                         recruitment = list(type = "pulse", years = 2012, mult = 8)),
             decline = treatment("recruitment -10%/yr",
                         recruitment = list(type = "trend", from = 2005, rate = -0.1))),
        nsim = 10, oem = oem_settings(one_gear_fishery(1), tail = 30), iters = 1:5)
plot_treatment(rt)
plot_length_indicators(length_indicators(rt, "Lutjanus bohar"))
```

### Main functions

| Step | Functions |
|---|---|
| Fishery | `fishery()`, `effort_trend()`, `reef_fishery()`, `load_fishery()`, `plot_effort()` |
| Species | `species()`, `selectivity()`, `species_library()`, `load_species()`, `save_species()` |
| Operating model | `build_om()`, `build_oms()`, `catch_fleet()`, `plot_om()`, `plot_fleets()` |
| Sampling | `oem_settings()`, `survey_settings()`, `season_design()`, `reef_oem()`, `sample_data()`, `lfd_annual()`, `plot_seasons()` (CPUE with AR1 residuals, FLRef two-stage length samples, biological samples with random or length-stratified ageing) |
| Colours | `floss.r4col()` (r4ss palette, as in JABBA) |
| Observations | `plot_availability()`, `plot_obs()`, `plot_lfd()`, `plot_bio()`, `alk_obs()`, `caa_obs()`, `plot_caa()` – observations only, nothing fitted |
| Treatments | `treatment()`, `run_treatment()`, `plot_treatment()`, `one_gear_fishery()`, `one_gear_species()` – baseline and treatment share all random numbers |
| Dynamics | `plot_eq()` (equilibrium curves with the realised path), `plot_srr()`, `plot_eq_selectivity()` (growth overfishing), `plot_ssb_vb()` (spawners vs what each gear sees) |
| Size data and flags | `length_indicators()`, `plot_length_indicators()`, `plot_lfd(cohorts = TRUE)`, `plot_indicators()` (dashboard), planted artefacts (`cpue_q_step`, `catch_bias_from`, `bio_age_bias`, survey `q_step`) with `planted_artefacts()` |
| Basic fits (optional) | `fit_growth()`, `fit_lw()`, `fit_maturity()`, `catch_curve()`, `catch_curve_length()`, `M_empirical()`, `cpue_check()`, `plot_fit()` |
| Truth (trainers) | `status_table()`, `plot_truth()`, `plot_kobe()`, `export_truth()` |
| Save | `export_data()`, `read_dataset()`, `simulate_fishery()` (all in one call) |

## Roadmap

| Version | Content |
|---|---|
| **0.1** | Specs (fishery, species) with YAML library; OM engine reproducing the multi-gear reef example; shared effort and effort creep |
| **0.2** | Sampling: catch, effort, CPUE (AR1 residuals), survey, two-stage length samples (FLRef `lfd.sim()`), biological samples (random or length-stratified ageing); observed ALK and catch-at-age; export with answer key; user guide |
| **0.3** | Treatments with common random numbers; dynamics plots; length indicators, cohort lines, indicator dashboard; plantable artefacts with hidden answer key; optional basic fits; guide and teaching vignette |
| **0.4** | Shiny app `run_app()` in the FLicc layout: Guide, Fishery, Species, OM, Treatments (predict, reveal), Sampling, Export, True status; info tooltips (`floss_tip()`), Hide truth switch |
| **0.4.1** | Seasonal sampling by default (`season_design()`: balanced, fixed imbalance, shift, opportunistic, per fleet), `lfd_annual()`, `plot_seasons()`, Seasons tab in the app; r4ss/JABBA palette `floss.r4col()` in all plots |
| **0.4.2** | Gradual sampling shifts, change in sampling intensity (e.g. funding cut), randomly missing sampling years |
| 0.4.x | Exporters for FLicc and JABBA3 |
| 0.5 | Trainee data viewer for exported datasets |
| 0.6 | Project save/load, more survey designs, seasonal effort |
| later | sex-structured and spatial stocks, link to `mse` via `sampling.lfd.oem()` |

## Acknowledgement

The progression of the teaching material (simple population dynamics, then data, then models) is inspired by Haddon, M. (2021) *Using R for Modelling and Quantitative Methods in Fisheries*, Chapman & Hall/CRC.

## License

GPL (>= 3)
