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

