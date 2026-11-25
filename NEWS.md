# FLOSS 0.6.1

* `amat_from_lm()`: age at 50% maturity from the length at maturity, inverting the growth curve. `species()` computes
  `tm` from `lm` when `tm` is missing or `NA`, so that the operating model (maturity at age) matures fish at `lm`.
* Shiny app: the age at 50% maturity is computed from the length at maturity, Linf, k and t0 and updates live
  (read-only), so maturity at age and at length always agree.
* Reef species (FishLife life histories; FishLife gives no t0):
  - t0 = -0.1 for the three teleosts (was -1, which put age-0 fish at 15-24 cm mid-year, at the peak of the trap
    selectivity, so age 0 was nearly fully selected). Blacktip reef shark keeps t0 = -1 (born large).
  - Age at maturity from FishLife's Lm for all species (FishLife's Lm and tm disagree): rabbitfish 1.56 (was 2.22),
    emperor 3.85 (4.85), snapper 3.01 (3.07), shark 5.83 (8.46).
  - Rabbitfish relq retuned (traps 0.58, gillnet 0.18, linefish 0.02; was 1.39, 0.42, 0.05) so that it stays just
    above Fmsy (F/Fmsy about 1.2 in 2025) as before; emperor and snapper trajectories are nearly unchanged.
  - Reproduce 0.6.0 OMs with the species files of FLOSS 0.6.0.
* `fishery(Escale)`: overall effort multiplier (default 1), applied to every fleet and year; `reef_3fleets` and
  `reef_fishery()` use 1.1, so all reef species are fished about 10% harder (2025 F/Fmsy: rabbitfish 1.29, emperor
  1.24, snapper 0.84). App: "Overall effort level" next to the effort scale between iterations.

# FLOSS 0.6.0

* Growth morphs (growth-type groups) in the operating model: `build_om(gtg = gtg_settings())` tracks morphs that
  differ in Linf through each cohort's fleet F history (length-selective mortality, Rosa Lee effect); the FLStock
  stays standard. Length and age samples come from the morphs in each fleet's catch. `gtg_settings(n, cv)` takes the
  total CV of Linf (within-morph CV = cv / sqrt(n), after M. Taylor's FLR platoon code). `gtg_lfq()` gives the true
  catch at length; `gtg_brp()` the morph-equilibrium reference points from steepness and R0. Plots:
  `plot_gtg_linf()` (design M. Taylor), `plot_gtg_growth()`, `plot_gtg_selectivity()`, `plot_gtg_lenage()`,
  `plot_gtg_lfq()`, `plot_gtg_spr()`. Guide: `floss_guide("growth-morphs")`.
* `fishery(Escale_cv, management)`: a lognormal effort scale per iteration (common to all fleets and species) and
  management changes (effort steps by fleet and year).
* New defaults (reproduce earlier OMs with FLOSS 0.5.8):
  - `reef_3fleets`: `Escale_cv = 0.25`; management change cutting gillnet and linefish effort to 50% and 55% in two
    steps (1998, 2002).
  - Rabbitfish relq traps 1.39, gillnet 0.42. Emperor relq gillnet 1.2; traps peak at 22 cm (sd_left 6, sd_right 15);
    linefish L50 38, L95 48 cm. Snapper relq gillnet 0.7, linefish 0.9; traps peak at 22 cm (sd_left 6).
    Recruitment variability (sigmaR) unchanged.
* `as_flr()` passes maturity L50 = L(tm) and L95 = L(tm + 1) as used by the OM.
* Guides render equations as MathML (no MathJax needed).

