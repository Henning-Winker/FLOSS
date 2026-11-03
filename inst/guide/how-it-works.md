## 1. What FLOSS does

FLOSS simulates fisheries whose truth is known and the data they produce.
It is built for teaching: trainees **predict** how a stock responds to a
change, see the response in the true dynamics and in the data, and learn to
recognise **warning flags**, patterns in data that real dynamics cannot
produce.

A FLOSS run has five layers:

1. **Fishery**: fleets, years and nominal effort, shared by all species.
2. **Species**: life history, target MSY, selectivity and catchability by
   fleet.
3. **Operating model (OM)**: the true age-structured dynamics.
4. **Sampling**: the data an observer would collect (catch, effort, CPUE,
   survey, length samples, biological samples).
5. **Truth and export**: one dataset per iteration, with a separate answer
   key.

A **treatment** changes one thing in the OM (F, recruitment, selectivity,
creep, M) and is run against an untreated baseline with the same random
numbers, so every difference comes from the treatment.

FLOSS describes dynamics and data. Fitting assessment models to the data
(e.g. with FLicc or JABBA) follows afterwards, using the exported data and
the answer key to judge the results.

## 2. Fishery and effort

Each fleet has a nominal effort trend, `E[f, y]`, with lognormal process
noise that is **shared by all species**: when a fleet fishes more in a year,
it fishes more on every species it catches. Effort creep raises the
efficiency of a fleet by a constant percentage per year; it increases
fishing mortality but is **not visible** in the nominal effort index.

    F[s, f, y] = relq[s, f] x Fmsy[s] x efficiency[f, y] x E[f, y]

`relq` is the relative catchability of species `s` by fleet `f`, set in the
species. With one fleet and `one_gear_species()`, effort equals F/Fmsy.

## 3. Species

A species holds von Bertalanffy growth (`linf`, `k`, `t0`), the
length-weight relationship (`a`, `b`; cm and g), maturity (`lm`, `tm`),
maximum age (also the plus group), natural mortality `M`, steepness `s`,
recruitment variability `sigR`, a target MSY, and selectivity and relative
catchability by fleet. Selectivity is defined on length (logistic, normal or
double normal) and converted to age through growth.

The four example species (*Siganus sutor*, *Lethrinus mahsena*, *Lutjanus
bohar*, *Carcharhinus melanopterus*) are read-only; copies can be saved to
the user's species library.

## 4. Operating model

The OM is an age-structured model built with FLR:

* the equilibrium is calibrated so that MSY equals the target (the unfished
  biomass is rescaled);
* fleet F-at-age is apical F times fleet selectivity; total F is the sum
  over fleets;
* the stock is projected with Beverton-Holt recruitment and lognormal
  deviations (optionally autocorrelated);
* reference points (Fmsy, Bmsy, MSY, B0) are recomputed for the realised
  fleet mix of the last years.

The **equilibrium curves** (yield, SSB and recruitment against F) describe
where the stock would settle under constant fishing. The realised years
usually lie off these curves: during a fishing-down phase, catches exceed
the equilibrium yield because biomass is being removed; after a reduction
in F the stock takes years to rebuild.

## 5. Treatments

| Treatment | What changes | Typical lesson |
|---|---|---|
| F step (`effort`) | effort of a fleet x `mult` from a year | lags: biomass and length structure respond over years |
| Fleet leaves (`effort`, `mult = 0`) | one fleet stops | other fleets' indicators and the fleet mix change |
| Recruitment pulse | recruitment x `mult` in some years | a cohort moves through the length data; catches rise later |
| Recruitment trend or regime | recruitment declines or shifts | mean length rises while the stock declines |
| Selectivity change | a fleet's gear changes from a year | length data and catch curves change without a change in F |
| Effort creep | efficiency rises from a year | CPUE stays high while the stock declines |
| M change | natural mortality x `mult` | productivity and reference points change |

Baseline and treated runs share every random number (effort noise,
recruitment deviations, sampling). Before the treatment year they are
identical.

## 6. Sampling

| Data | How it is generated | Key settings |
|---|---|---|
| Catch by fleet | true catch x (1 - under-reporting) x lognormal error | CV, bias, bias from year |
| Effort index | nominal effort x lognormal error | CV |
| CPUE by fleet | vulnerable biomass^beta x efficiency (if creep is carried) x lognormal AR1 error | CV, rho, beta, catchability step |
| Survey index | biomass at survey time x survey selectivity, AR1 error | timing, CV, rho, selectivity (relative to Linf), catchability step |
| Length samples | two-stage sampling by season: clusters drawn from the season's catch-at-age, then fish measured through the gear's length key at the middle of the season | clusters (`ess_age`), fish measured (`ess_len`), length bin, seasonal design |
| Biological samples | one row per fish and season: length (mm), age (with ageing error and bias), weight (g), maturity (0/1) | fish per year, ageing scheme (random or length-stratified), CVs, seasonal design |

CPUE, length and biological samples cover the last years of each fleet
(e.g. 5, 10 and 15 years in the reef example); catch and effort cover all
years.

**Two-stage length sampling.** Fish in a trip or haul are alike, so a length
sample carries less information than its number of fish suggests. The
effective sample size lies between the number of clusters and the number
of fish measured, depending on growth and selectivity.

**Seasonal sampling.** The stock is annual, with fishing effort spread
evenly through the year; only the sampling is seasonal. By default the length and biological samples of a
year are spread evenly over four seasons (`season_design()`). Within a
season the catch-at-age is the catch taken between the start and end of
the season, and lengths are those at the middle of the season, so fish in
late-season samples are larger. Imbalanced designs:

| Design | Example | Effect on the pooled annual length data |
|---|---|---|
| Balanced | equal shares | none |
| Fixed imbalance | most samples in season 1 | mean length too low every year (a constant bias) |
| Shift | sampling moves to seasons 3 and 4 from 2018 | a step in mean length in 2018 that looks like a change in mortality |
| Opportunistic | shares vary at random from year to year | extra year-to-year noise in mean length |
| Missing season | share 0 in the monsoon | no information on that part of the year |
| Gradual shift | the shift spread over several years | a slow drift in mean length, easily read as a trend in the stock |
| Fewer samples | sample size x 0.3 from 2020 (funding cut) | noisier length data, more erratic indicators |
| Missing years | each year unsampled with some probability | gaps; cohorts harder to follow |

`lfd_annual(weight = "equal")` and `length_indicators(season_weight =
"equal")` give every sampled season the same weight, which corrects an
uneven spread but cannot recover a season that was not sampled.
`plot_seasons()` shows the sample sizes and the mean length by season.

**Autocorrelated CPUE residuals.** A fleet that keeps fishing the same
grounds does not sample the stock at random; its errors persist from year
to year (`rho` > 0), so each year adds less independent information.

## 7. Size data: what to expect

At equilibrium (constant recruitment and mortality, flat-topped
selectivity) a length composition has a fixed shape:

* the **ascending limb** is set by selectivity (`Lc`, the length of first
  capture);
* the **descending limb** is set by mortality relative to growth, `Z/k`;
* **mean length** above `Lc` follows Beverton and Holt (1956):

      Lmean = (Linf + (Z/k) Lc) / (1 + Z/k)    so    Z = k (Linf - Lmean) / (Lmean - Lc)

Heavier fishing gives a lower mean length, a shorter right-hand tail and
fewer large fish (lower `L95`, `Pmega`).

What breaks the pattern:

| Violation | What the length data show | Misread as | Clue |
|---|---|---|---|
| F has just changed | mean length adjusts over years | F still high (or low) | indicators still trending |
| Strong cohort | a mode enters small and moves right | first overfishing, then recovery | the mode moves at the growth rate; other indices show the year class |
| Recruitment declining | fewer small fish, mean length rises | recovery | juvenile index and CPUE fall |
| Recruitment rising | more small fish, mean length falls | overfishing | juvenile index rises |
| Dome-shaped gear | large fish missing whatever F | high F | a logistic gear or the survey still sees large fish |
| Gear or target change | step in length composition | sudden change in F | coincides with a fleet event |
| Low effective sample size | erratic year-to-year shapes | recruitment variation | consecutive years not consistent with growth |

Real modes move to the right by about one year's growth per year, fade at
rate Z, and appear in every gear that can catch them. Patterns that jump,
vanish or appear in one gear only are sampling noise or selectivity.

## 8. Warning flags

| Flag in the data | Likely cause, not real dynamics |
|---|---|
| CPUE flat or rising while survey, catch and mean length decline | effort creep or hyperstable CPUE |
| Catch well above MSY for years with no decline in any index | over-reported catch, or indices that do not track abundance |
| Mean length rises while F rises | dome-shaped gear or a gear change |
| Catch-curve Z below M | selectivity not flat on the descending limb, or a recruitment trend |
| A strong year class in an index that never moves through the length data | index noise or ageing error |
| A step in an index in one year | change in catchability (vessel, gear, electronics) |
| A step in mean length in the year the sampling season changed | seasonal sampling shift (fish are larger later in the year) |
| CPUE residuals running in long streaks | non-random sampling |
| CPUE of a gear that catches juveniles (or lets large fish escape) declines much less than other indices | the gear's vulnerable biomass is not the spawning stock (see `plot_ssb_vb()`) |

Plausibility checks: a stock cannot grow faster than its productivity
allows; cohorts must move through lengths at the growth rate; biomass does
not respond to F instantly; and independent indicators of the same stock
should broadly agree. When they do not, at least one of them is not
tracking the stock.

## 9. Truth and export

`export_data()` writes one folder (or workbook) per iteration with the
observed data and a README, plus a separate `answer_key/` folder with the
true time series, reference points, true status, true catch, true ages and
the list of **planted artefacts** (data problems and treatments set on
purpose). Remove the answer key before giving the data to trainees.
`export_truth()` writes the truth alone, for the trainer.
