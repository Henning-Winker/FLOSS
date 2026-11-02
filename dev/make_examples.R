# Builds the example species and fishery YAML files in inst/extdata from the
# multi-gear reef operating model (FLsim_multigearreef_revised).
# Run from the package root: Rscript dev/make_examples.R

devtools_load <- function() {
  if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".") else
    library(FLOSS)
}
devtools_load()

load("dev/fishlife_lhpars.rdata")   # lhpars: FishLife traits via SPMpriors

msy <- c(siganus_sutor = 750, lethrinus_mahsena = 1000, lutjanus_bohar = 600,
         carcharhinus_melanopterus = 75)
common <- c(siganus_sutor = "Whitespotted rabbitfish",
            lethrinus_mahsena = "Mahsena emperor",
            lutjanus_bohar = "Two-spot red snapper",
            carcharhinus_melanopterus = "Blacktip reef shark")
prod <- c(siganus_sutor = "High", lethrinus_mahsena = "Intermediate",
          lutjanus_bohar = "Low-intermediate",
          carcharhinus_melanopterus = "Very low")

sel <- list(
  siganus_sutor = list(
    traps    = selectivity("dnormal", Lpeak = 20, sd_left = 5, sd_right = 15),
    gillnet  = selectivity("normal", Lpeak = 32, sd = 8),
    linefish = selectivity("logistic", L50 = 28, L95 = 32)),
  lethrinus_mahsena = list(
    traps    = selectivity("dnormal", Lpeak = 27, sd_left = 8, sd_right = 0.7 * 28),
    gillnet  = selectivity("normal", Lpeak = 35, sd = 0.3 * 35),
    linefish = selectivity("logistic", L50 = 35, L95 = 50)),
  lutjanus_bohar = list(
    traps    = selectivity("dnormal", Lpeak = 27, sd_left = 10, sd_right = 0.5 * 27),
    gillnet  = selectivity("normal", Lpeak = 37, sd = 0.3 * 37),
    linefish = selectivity("logistic", L50 = 40, L95 = 55)),
  carcharhinus_melanopterus = list(
    traps    = selectivity("dnormal", Lpeak = 50, sd_left = 15, sd_right = 50),
    gillnet  = selectivity("dnormal", Lpeak = 70, sd_left = 20, sd_right = 80),
    linefish = selectivity("logistic", L50 = 70, L95 = 110)))

relq <- list(
  siganus_sutor             = c(traps = 1.00, gillnet = 0.30, linefish = 0.05),
  lethrinus_mahsena         = c(traps = 0.30, gillnet = 0.70, linefish = 1.00),
  lutjanus_bohar            = c(traps = 0.10, gillnet = 0.50, linefish = 1.00),
  carcharhinus_melanopterus = c(traps = 0.02, gillnet = 1.40, linefish = 0.40))

dir.create("inst/extdata/species", recursive = TRUE, showWarnings = FALSE)
for (i in seq_len(nrow(lhpars))) {
  lh <- as.list(lhpars[i, -1])
  id <- FLOSS:::species_id(lhpars$species[i])
  sp <- species(species = lhpars$species[i], common_name = common[[id]],
                productivity = prod[[id]], lh = lh, msy = msy[[id]],
                sel = sel[[id]], relq = relq[[id]])
  write_spec(sp, file.path("inst/extdata/species", paste0(id, ".yaml")))
}

dir.create("inst/extdata/fishery", recursive = TRUE, showWarnings = FALSE)
write_spec(reef_fishery(), "inst/extdata/fishery/reef_3fleets.yaml")
message("written: ", toString(list.files("inst/extdata", recursive = TRUE)))
