# FLOSS <img src="https://img.shields.io/badge/status-v0.2%20prototype-orange" align="right"/>

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
* **Engine first.** Everything runs from R scripts; the Shiny app (v0.3)
  only collects inputs and calls the same functions.

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

The full walk-through is the vignette *FLOSS user guide*
(`vignette("floss-02-user-guide", package = "FLOSS")`).

### Main functions

| Step | Functions |
|---|---|
| Fishery | `fishery()`, `effort_trend()`, `reef_fishery()`, `load_fishery()`, `plot_effort()` |
| Species | `species()`, `selectivity()`, `species_library()`, `load_species()`, `save_species()` |
| Operating model | `build_om()`, `build_oms()`, `catch_fleet()`, `plot_om()`, `plot_fleets()` |
| Sampling | `oem_settings()`, `survey_settings()`, `reef_oem()`, `sample_data()` (CPUE with AR1 residuals, FLRef two-stage length samples, biological samples with random or length-stratified ageing) |
| Observations | `plot_availability()`, `plot_obs()`, `plot_lfd()`, `plot_bio()`, `alk_obs()`, `caa_obs()`, `plot_caa()` – observations only, nothing fitted |
| Truth (trainers) | `status_table()`, `plot_truth()`, `plot_kobe()`, `export_truth()` |
| Save | `export_data()`, `read_dataset()`, `simulate_fishery()` (all in one call) |

## Roadmap

| Version | Content |
|---|---|
| **0.1** | Specs (fishery, species) with YAML library; OM engine reproducing the multi-gear reef example; shared effort and effort creep |
| **0.2** | Sampling: catch, effort, CPUE (AR1 residuals), survey, two-stage length samples (FLRef `lfd.sim()`), biological samples (random or length-stratified ageing); observed ALK and catch-at-age; export with answer key; user guide |
| 0.3 | Shiny app: Fishery - Species - OM - Sampling - Export, with a "hide truth" switch for live teaching |
| 0.4 | Trainee data viewer for exported datasets |
| 0.5 | Project save/load, more survey designs, seasonal effort |
| later | sex-structured and spatial stocks, link to `mse` via `sampling.lfd.oem()` |

## License

GPL (>= 3)
