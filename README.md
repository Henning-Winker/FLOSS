# FLOSS <img src="https://img.shields.io/badge/status-v0.1%20prototype-orange" align="right"/>

**FL Operating model & Sampling Simulator**

FLOSS builds multi-fleet, multi-species age-structured operating models and
(from v0.2) the data such fisheries produce: catch and effort by fleet, CPUE,
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

fy <- load_fishery("reef_3fleets")          # traps, gillnet, linefish 1969-2025
species_library()                           # example + user species

oms <- build_oms(species_library()$id, fy, nsim = 100)
oms$lutjanus_bohar
plot_om(oms)                                # B/Bmsy, F/Fmsy, C/MSY by species
plot_fleets(oms$lutjanus_bohar)             # apical F by fleet
catch_fleet(oms$lutjanus_bohar)             # true catch by fleet

# a new species, saved to your library
sp <- load_species("lutjanus_bohar")
sp <- species("Lutjanus gibbus", common_name = "Humpback red snapper",
              lh = modifyList(sp$lh, list(linf = 45, k = 0.3, tmax = 18)),
              msy = 400, sel = sp$sel, relq = c(traps = 0.4, gillnet = 0.6,
                                                linefish = 0.8))
save_species(sp)

# effort creep of 2% per year in the linefish fleet
fy2 <- reef_fishery(creep = c(traps = 0, gillnet = 0, linefish = 0.02))
plot_effort(fy2, efficiency = TRUE)
```

## Roadmap

| Version | Content |
|---|---|
| **0.1** | Specs (fishery, species) with YAML library; OM engine reproducing the multi-gear reef example; shared effort and effort creep |
| 0.2 | Sampling (OEM): catch, effort, CPUE, survey, LFD (via `FLRef::sample_lfd_gears()`), biological samples; export tables + answer key |
| 0.3 | Shiny app: Fishery - Species - OM - Sampling - Export, trainer/trainee mode |
| 0.4 | Project save/load, more survey designs, seasonal effort |
| later | sex-structured and spatial stocks, link to `mse` via `sampling.lfd.oem()` |

## License

GPL (>= 3)
