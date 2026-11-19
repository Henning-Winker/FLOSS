# FLOSS app guide

<!-- App-only copy in plain language. The website article inst/guide/how-it-works.md is unchanged. Sections are found by their heading, so headings must stay as they are. -->

## 1. Fishery and effort
Each fleet follows a nominal effort trend: effort starts at a low level in the first year, rises to a peak, moves to a later level, and then remains there. Around this trend, effort varies randomly from year to year, and this variation is shared by all species a fleet catches, so when a fleet fishes more in a given year, it fishes more on every species in its catch. Effort creep raises a fleet's fishing efficiency by a constant percentage each year, reflecting improvements in electronics, gear, or skill; it increases fishing mortality but does not appear in the nominal effort that fishers report.

The fishing mortality a fleet exerts on a species in a given year is the product of four terms: the relative catchability of that species by the fleet, the species' fishing mortality at MSY (Fmsy), the fleet's fishing efficiency in that year, and the fleet's effort in that year. Relative catchability, set in the Species step, determines how strongly a fleet fishes each species; with relative catchability and effort both equal to 1, and no effort creep, a fleet fishing alone exploits the species at approximately Fmsy. In the one-gear fishery, catchability is set so that effort equals the ratio of fishing mortality to Fmsy (F/Fmsy) directly.

The New fishery option starts with no gears: users add gears by name and set each gear's effort. Every species then receives a default S-shaped selectivity for each new gear, centred on its length at maturity, which can be adjusted in the Species step.

## 2. Species
Each species is defined by its von Bertalanffy growth parameters (asymptotic length, Linf; growth rate, k; and theoretical age at zero length, t0); its length-weight relationship (a and b, with length in cm and weight in g); its length and age at 50% maturity; its maximum age, which also forms the plus group; natural mortality (M); the steepness of its stock-recruitment relationship; the variability of its recruitment; a target MSY; and its selectivity and relative catchability by fleet. Selectivity is defined on length, as an S-shaped (logistic), symmetric dome-shaped (normal), or asymmetric dome-shaped (double normal) curve, and is converted to age through the growth curve.

The four example species (Siganus sutor, Lethrinus mahsena, Lutjanus bohar, and Carcharhinus melanopterus) are read-only, although copies can be saved to the user's own species library and edited there.

## 3. Operating model
The operating model represents the true, age-structured dynamics of each stock and is built with the FLR fisheries modelling framework:
- each stock is scaled so its equilibrium MSY equals the target MSY, by rescaling its unfished biomass;
- each fleet's fishing mortality at age equals its apical (highest) fishing mortality multiplied by its selectivity at age, and total fishing mortality is the sum across all fleets;
- the stock is projected year by year with Beverton-Holt recruitment and random lognormal deviations, which can optionally be autocorrelated, so good and poor recruitment years cluster together;
- reference points (Fmsy, the spawning biomass at MSY, Bmsy, MSY, and unfished biomass, B0) are recalculated for the realised fleet mix of the final years.

The equilibrium curves (yield, spawning biomass, and recruitment plotted against fishing mortality) describe where the stock would settle if fished at a constant level for a long time. The simulated years usually lie off these curves: during a fishing-down phase, catches exceed the equilibrium yield because accumulated biomass is being removed, while after a reduction in fishing mortality, the stock takes years to rebuild.

## 4. Sampling
Steps 1 to 3 form the first part of the app: they build the true stock, which in a real fishery could never be observed directly. Step 4 begins the second part, in which FLOSS generates the data a monitoring programme would collect from that stock. For each fleet and year, it draws reported catch and effort, catch per unit effort (CPUE), an optional fishery-independent survey index, length samples, and biological samples recording the length, age, weight, and maturity of individual fish. Like real data, these observations are deliberately imperfect: each carries the observation error, bias, and sampling design set in the Sampling step, so users control how far the data depart from the truth. Each iteration of the operating model can yield one dataset, and the table below summarises how each data type is generated and which settings control it.

| Data | How it is generated | Key settings |
|---|---|---|
| Catch by fleet | the true catch, reduced by any under-reporting, with random observation error | observation error (CV), percentage under-reported, and the year under-reporting begins |
| Effort by fleet | the nominal effort, reported in fishing days, with random observation error; effort creep is never included | observation error (CV) and fishing days at a relative effort of 1 |
| CPUE by fleet | follows the biomass vulnerable to the fleet, optionally falling more slowly than that biomass (hyperstability), rising with the fleet's efficiency when effort creep is carried into the index, and carrying random errors that can persist from year to year | observation error (CV), autocorrelation, proportionality, and catchability step |
| Survey index (optional) | follows the biomass selected by the survey at the time of year the survey takes place, with random errors that can persist from year to year | survey timing, observation error (CV), autocorrelation, survey selectivity (set relative to Linf), and catchability step |
| Length samples | drawn in two stages each season: clusters (trips or hauls) are sampled from the season's catch, and fish are then measured from those clusters at their mid-season length | clusters per year, fish measured per year, length bin, and seasonal design |
| Biological samples | one record per fish, holding its length, its age (with optional ageing error and bias), its weight, and whether it is mature | fish per year, ageing scheme (a random share or length-stratified), ageing and weighing error (CV), and seasonal design |

Lengths are entered in centimetres in the app, whereas the exported length and biological samples record length in millimetres and weight in grams.

Catch and effort are reported for every year, whereas CPUE, length samples, and biological samples cover only the most recent years of each fleet (for example, 5, 10, and 15 years for the three fleets in the reef example), as is typical of monitoring programmes that began after the fishery developed.

**Two-stage length sampling.** Fish caught on the same trip or in the same haul tend to resemble one another, so a length sample carries less information than its number of fish suggests. Its effective sample size lies between the number of clusters and the number of fish measured, depending on growth and selectivity.

**Seasonal sampling.** The stock itself is modelled annually, with fishing effort spread evenly through the year; only the sampling is seasonal. By default, the length and biological samples of each year are spread evenly over four seasons. The fish sampled in a season come from the catch taken during that season and are measured at their mid-season length, so fish in late-season samples are larger. Imbalanced designs, set in the Seasons tab of the Sampling step, affect the pooled annual length data as follows:

| Design | Example | Effect on the pooled annual length data | How to set it in the Seasons tab |
|---|---|---|---|
| Balanced | equal shares in every season | none | Design: Balanced |
| Fixed imbalance | most samples taken in season 1 | mean length too low in every year (a constant bias) | Design: Fixed imbalance, with the share sliders |
| Shift | sampling moves to seasons 3 and 4 from 2018 | a step in mean length in 2018 that resembles a change in mortality | Design: Shift, with Shift from year and the second set of share sliders |
| Opportunistic | shares vary at random from year to year | additional year-to-year noise in mean length | Design: Opportunistic, with the average share sliders and How regular |
| Missing season | no samples during the monsoon | no information on that part of the year | Fixed imbalance or Shift, with the share of that season set to 0 |
| Gradual shift | the shift spread over several years | a slow drift in mean length, easily mistaken for a trend in the stock | Design: Shift, with Gradually over set above 0 |
| Fewer samples | sample size reduced to 30% from 2020 (for example, after a funding cut) | noisier length data and more erratic indicators | Sample size changes from a given year, with the sample size multiplier |
| Missing years | each year left unsampled with a set probability | gaps in the series, making cohorts harder to follow | Chance that a year is not sampled |

Weighting every sampled season equally when pooling a year's samples corrects an uneven spread, but cannot recover a season that was not sampled. The Seasons tab plots sample sizes and mean length by season, with the annual mean shown both ways.

**Autocorrelated CPUE errors.** A fleet that keeps returning to the same fishing grounds does not sample the stock at random, so its errors persist from year to year and each additional year adds less independent information. This persistence is set with the autocorrelation slider in the CPUE tab of the Sampling step.

### Reading the length data
*Where this applies in the app: Sampling, Length frequencies tab (length samples and indicators) and Seasons tab; and Treatments (step 5), Reveal tab ("What the length data show").*

Under equilibrium conditions, meaning constant recruitment and mortality and a gear that retains all fish above a certain size, a length distribution takes a fixed shape:
- the rising left-hand side is set by selectivity, through the length at first capture (Lc, taken as the most common length in the sample);
- the falling right-hand side is set by total mortality (Z) relative to the growth rate (k);
- the mean length of fish at or above Lc follows Beverton and Holt (1956): mean length = (Linf + (Z/k) × Lc) / (1 + Z/k), which can be rearranged to estimate total mortality as Z = k × (Linf − mean length) / (mean length − Lc).

Heavier fishing therefore produces a lower mean length, a shorter right-hand tail, and fewer large fish, shown by a lower 95th-percentile length (L95) and a smaller proportion of "mega-spawners" (Pmega, the share of fish more than 10% longer than the optimum length at which a cohort's biomass peaks).

The pattern breaks down whenever these conditions do not hold:

| Situation | What the length data show | Commonly misread as | Clue |
|---|---|---|---|
| Fishing mortality has just changed | mean length adjusts gradually over several years | fishing mortality still high (or still low) | indicators are still trending |
| Strong cohort | a mode enters at small lengths and moves to the right | first overfishing, then recovery | the mode moves at the growth rate, and other indices show the same year class |
| Declining recruitment | fewer small fish, so mean length rises | recovery | the juvenile index and CPUE fall |
| Rising recruitment | more small fish, so mean length falls | overfishing | the juvenile index rises |
| Dome-shaped gear | large fish missing regardless of fishing mortality | high fishing mortality | an S-shaped gear or the survey still catches large fish |
| Change of gear or target | a step in the length composition | a sudden change in fishing mortality | the step coincides with a fleet event |
| Low effective sample size | erratic shapes from year to year | recruitment variation | consecutive years are inconsistent with growth |

Real modes move to the right by about one year's growth each year, fade at the rate of total mortality, and appear in every gear able to catch them; patterns that jump, vanish, or appear in a single gear reflect sampling noise or selectivity.

## 5. Treatments
Treatments are optional. Each one changes a single element of the operating model from a chosen year and is compared with the baseline built in step 3, while its length data are sampled with the settings chosen in step 4, so the Sampling settings should be set before a treatment is run. Up to four treatments can be saved and compared within the app, although their data are not exported.

| Treatment | What changes | Typical lesson |
|---|---|---|
| Fishing effort step | a fleet's effort is multiplied by a chosen factor from a given year | responses lag: biomass and length structure change over several years |
| Fleet leaves | one fleet's effort falls to zero | indicators from the remaining fleets, and the fleet mix, change |
| Recruitment pulse | recruitment is multiplied in one or more years | a strong or weak cohort moves through the length data, and catches rise or fall later |
| Recruitment trend or regime shift | recruitment declines steadily, or shifts to a new level | mean length can rise while the stock declines |
| Selectivity change | a fleet's gear changes from a given year | length data and catch curves change without any change in fishing mortality |
| Effort creep | a fleet's efficiency rises each year from a given year | CPUE stays high while the stock declines |
| Natural mortality change | natural mortality is multiplied by a chosen factor | productivity and reference points change |

Baseline and treated runs share every random number (effort noise, recruitment deviations, and sampling), so the two are identical until the treatment year, and every later difference results from the treatment.

## Warning flags
The table below links patterns users may find in the simulated data to their most likely cause. Each of these problems can be planted deliberately in the Fishery and Sampling steps, and those planted in a run are listed in the Planted warning flags box on the Overview tab of the Sampling step (shown when Hide truth is unticked) and in the answer key of each export.

| Flag in the data | Likely cause, not real dynamics |
|---|---|
| Catch per unit effort (CPUE) stays flat or rises while the survey index, catch, and mean length decline | effort creep, or hyperstable CPUE (an index that falls more slowly than the biomass it tracks) |
| Catch remains well above maximum sustainable yield (MSY) for years, with no decline in any index | over-reported catch, or indices that do not track abundance |
| Mean length rises while fishing mortality (F) rises | a dome-shaped gear, or a change of gear |
| Total mortality (Z) estimated from a catch curve, which measures how quickly numbers decline with age or length in the catch, falls below natural mortality (M) | selectivity that is not flat across the larger fish, or a trend in recruitment |
| A strong year class appears in an index but never moves through the length data | noise in the index, or ageing error |
| An index steps up or down in a single year | a change in catchability (for example, a new vessel, gear, or electronics) |
| Mean length steps up or down in the year the sampling season changed | a shift in seasonal sampling, since fish are larger later in the year |
| CPUE residuals (the differences between the index and its underlying trend) run in long streaks | non-random sampling, such as a fleet returning to the same fishing grounds |
| CPUE from a gear that catches juveniles, or lets large fish escape, declines much less than other indices | the gear's vulnerable biomass (the part of the stock it can catch) is not the spawning stock; the two can be compared in the Spawners vs vulnerable tab of the Operating model step |

**Plausibility checks.** A stock cannot grow faster than its productivity allows; cohorts must move through the length data at the growth rate; biomass does not respond to a change in fishing mortality instantly; and independent indicators of the same stock should broadly agree. When they do not, at least one indicator is not tracking the stock.

## 6. Export
Downloading the simulated datasets allows them to be analysed outside the app, for example by fitting stock assessment models in FLicc or JABBA, and keeps a fixed set of exercises that trainees can work through and trainers can later check against the truth. To download, users open the Datasets tab of the Export step, choose the species, the datasets (iterations), the information to give trainees, and the file format, and then press Download datasets (.zip); the Settings tab of the same step saves the fishery, species, and sampling settings, with the seed, so the same data can be re-created.

Each download contains one folder of comma-separated (csv) files, or one Excel workbook, per dataset, holding the observed data and a README file that describes them. A separate answer key folder holds the true time series, the reference points, the true stock status, the true catch, the true age of every sampled fish, and a list of planted artefacts (the data problems and treatments set deliberately). The answer key can be included only when Hide truth is unticked and its tick box in the Export step is ticked, and it should be left out of any download intended for trainees.

**Naive aggregation (teaching).** Because effort is reported in fishing days, catch, effort, and nominal catch per unit effort (CPUE, in kilograms per fishing day) can be summed across fleets, or across both fleets and species, and length samples can be pooled across fleets, either as measured or raised to each fleet's catch. These aggregates are shown in the Catch, Effort, CPUE, and Length frequencies tabs of the Sampling step, and ticking the option for aggregated tables in the Export step adds them to the download. Biological samples can likewise be grouped by fleet, across all fleets, or by season in the Biological samples tab of the Sampling step.

| Aggregate | What it ignores | What it can produce |
|---|---|---|
| Effort summed across fleets | a fishing day differs between gears | effort that rises while fishing mortality falls, as the fleet mix changes |
| Nominal CPUE across all fleets | differences in catchability and selectivity between fleets | trends driven by changes in the fleet mix |
| Nominal CPUE across all species | species composition | a stable index while one species collapses |
| Length samples pooled as measured | each fleet's share of the catch | a length distribution dominated by the most heavily sampled fleet |

**Files for R users.** The R data (.rds) format stores the data as objects of the FLR fisheries modelling framework, holding catch, effort, and length frequencies by year, abundance indices with their standard errors on the log scale, and the species' life-history parameters, ready to load into the JABBA and FLicc assessment packages.

## 7. True status
Trainers can download the truth on its own, without any observed data, from the Downloads tab of the True status step: the true time series, the reference points, and the true stock status for every iteration of the operating model.
