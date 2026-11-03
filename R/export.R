# Wrappers, export and import of datasets.

#' Simulate a fishery and its data in one call
#'
#' Wrapper around [build_oms()] and [sample_data()].
#'
#' @param species species ids from the library or [species()] objects.
#' @param fishery a [fishery()].
#' @param oem [oem_settings()]; default `oem_settings(fishery)`.
#' @param nsim iterations of the operating models.
#' @param iters iterations to sample as datasets.
#' @param seed base seed for sampling.
#' @param ... passed to [build_oms()] (e.g. `sigmaQ`, `ar1`).
#' @return a list with `oms` (`floss_oms`) and `data` (`floss_data`).
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_fishery("siganus_sutor", reef_fishery(), nsim = 5, iters = 1)
#' sim$data
#' }
simulate_fishery <- function(species = species_library()$id,
                             fishery = load_fishery(), oem = NULL, nsim = 20,
                             iters = 1:min(10, nsim), seed = 123, ...) {
  if (is.null(oem)) oem <- oem_settings(fishery)
  oms <- build_oms(species, fishery, nsim = nsim, ...)
  list(oms = oms, data = sample_data(oms, oem, iters = iters, seed = seed))
}

#' Sampling settings of the reef example
#'
#' Length, CPUE and biological samples cover the most recent 5 years for
#' traps, 10 for gillnet and 15 for linefish; linefish CPUE residuals are
#' more autocorrelated (clustered fishing).
#' @param ... further arguments to [oem_settings()].
#' @return a `floss_oem`.
#' @export
reef_oem <- function(...) {
  args <- utils::modifyList(list(
    fishery = load_fishery("reef_3fleets"),
    tail = c(traps = 5, gillnet = 10, linefish = 15),
    cpue_rho = c(traps = 0.3, gillnet = 0.3, linefish = 0.6),
    effort_days = c(traps = 4000, gillnet = 1500, linefish = 3000),
    lfd_ess_age = c(traps = 30, gillnet = 30, linefish = 50),
    lfd_ess_len = c(traps = 300, gillnet = 300, linefish = 500),
    bio_n = c(traps = 100, gillnet = 100, linefish = 200)), list(...))
  do.call(oem_settings, args)
}

#' Export datasets for trainees
#'
#' Writes one folder (or workbook) per iteration with the observed data, the
#' settings needed to re-create them, and, separately, the answer key:
#'
#' ```
#' <dir>/README.txt, settings.yaml
#' <dir>/iter_01/catch.csv effort.csv cpue.csv survey.csv lfd.csv bio.csv [lhpars.csv]
#' <dir>/answer_key/README_answer_key.txt truth.csv status.csv refpts.csv
#'   catch_true.csv bio_true_age.csv planted_artefacts.csv
#' ```
#'
#' @param data a `floss_data`.
#' @param dir output folder (created).
#' @param give what goes into the trainee folders: `"data"` or `"data_lh"`
#'   (adds the life-history parameters).
#' @param format any of `"csv"`, `"xlsx"` (needs 'writexl'; one workbook per
#'   iteration) and `"rds"` (one `floss_flr.rds` per iteration: a list by
#'   species of [as_flr()] objects, i.e. FLQuants of catch, effort and length
#'   frequencies, FLIndices and the life-history FLPar, ready for JABBA3 and
#'   FLicc).
#' @param answer_key write the `answer_key` folder.
#' @param aggregate also write the aggregated tables (for teaching):
#'   `catch_aggregate` (by fleet, all fleets, all species), `effort_aggregate`
#'   (by fleet, all fleets), `cpue_nominal` (catch / effort at the three
#'   levels) and `lfd_aggregate` (by fleet and all fleets, as sampled and
#'   raised to the catch); see [catch_table()] and [lfd_aggregate()].
#' @param zip also write `<dir>.zip`.
#' @return the folder, invisibly.
#' @export
export_data <- function(data, dir, give = c("data_lh", "data"),
                        format = "csv", answer_key = TRUE, aggregate = FALSE,
                        zip = FALSE) {
  give <- match.arg(give)
  stopifnot(inherits(data, "floss_data"))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  its <- data$settings$iters
  streams <- c("catch", "effort", "cpue", "survey", "lfd", "bio")
  for (i in its) {
    tabs <- lapply(streams, function(s) {
      d <- data[[s]]
      if (is.null(d)) return(NULL)
      d <- d[d$iter == i, setdiff(names(d), c("iter", "catch_true",
                                              "effort_true", "index_true",
                                              "age_true")), drop = FALSE]
      rownames(d) <- NULL
      d
    })
    names(tabs) <- streams
    if (give == "data_lh") tabs$lhpars <- data$lhpars
    if (aggregate) tabs <- c(tabs, aggregate_tables(data, i))
    tabs <- Filter(Negate(is.null), tabs)
    sub <- file.path(dir, sprintf("iter_%02d", i))
    if ("csv" %in% format) {
      dir.create(sub, showWarnings = FALSE)
      for (s in names(tabs))
        utils::write.csv(tabs[[s]], file.path(sub, paste0(s, ".csv")),
                         row.names = FALSE)
    }
    if ("xlsx" %in% format) {
      if (!requireNamespace("writexl", quietly = TRUE))
        stop("format = 'xlsx' needs the 'writexl' package")
      writexl::write_xlsx(tabs, paste0(sub, ".xlsx"))
    }
    if ("rds" %in% format) {
      dir.create(sub, showWarnings = FALSE)
      flr <- as_flr(data, iter = i, all = TRUE)
      if (give == "data") flr <- lapply(flr, function(x) { x$lhpar <- NULL; x })
      saveRDS(flr, file.path(sub, "floss_flr.rds"))
    }
  }
  if (answer_key) {
    ak <- file.path(dir, "answer_key")
    dir.create(ak, showWarnings = FALSE)
    utils::write.csv(truth_long(data), file.path(ak, "truth.csv"), row.names = FALSE)
    utils::write.csv(status_table(data), file.path(ak, "status.csv"), row.names = FALSE)
    arte <- data$settings$artefacts
    utils::write.csv(arte, file.path(ak, "planted_artefacts.csv"), row.names = FALSE)
    writeLines(c("ANSWER KEY - do not give to trainees",
                 "",
                 "truth.csv              true SSB, recruitment, Fbar, catch, B/Bmsy, F/Fmsy by year",
                 "status.csv             true status in the last year (Kobe)",
                 "refpts.csv             reference points",
                 "catch_true.csv         true catch by fleet",
                 "bio_true_age.csv       true age of every sampled fish",
                 "planted_artefacts.csv  data problems and treatments planted on purpose",
                 "",
                 "Planted artefacts:",
                 if (nrow(arte)) paste0("  ", arte$artefact, " (", arte$where, "): ", arte$detail)
                 else "  none"),
               file.path(ak, "README_answer_key.txt"))
    utils::write.csv(data$refpts, file.path(ak, "refpts.csv"), row.names = FALSE)
    utils::write.csv(data$bio[, c("species", "iter", "year", "source",
                                  "fish_id", "age_true")],
                     file.path(ak, "bio_true_age.csv"), row.names = FALSE)
    obs_true <- function(d, key) if (!is.null(d)) d[, key, drop = FALSE]
    utils::write.csv(obs_true(data$catch, c("species", "iter", "year", "fleet",
                                            "catch_true")),
                     file.path(ak, "catch_true.csv"), row.names = FALSE)
  }
  st <- data$settings
  yaml::write_yaml(list(
    floss_version = st$floss_version, iters = st$iters, seed = st$seed,
    nsim = st$nsim, fishery = unclass_deep(st$fishery),
    species = lapply(st$species, unclass_deep), oem = unclass_deep(st$oem)),
    file.path(dir, "settings.yaml"))
  writeLines(readme_text(data, give, answer_key), file.path(dir, "README.txt"))
  if (zip) {
    old <- setwd(dirname(normalizePath(dir)))
    on.exit(setwd(old))
    utils::zip(paste0(basename(dir), ".zip"), basename(dir), flags = "-rq")
  }
  invisible(dir)
}

#' Read an exported dataset
#'
#' Reads the csv files of one exported iteration back into a list of data
#' frames, so that trainees (or the plotting functions) can use them.
#' @param dir the iteration folder (e.g. `".../iter_01"`).
#' @return a list of data frames (`catch`, `effort`, `cpue`, `survey`, `lfd`,
#'   `bio`, `lhpars` where present), each with an `iter` column set to 1.
#' @export
read_dataset <- function(dir) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  out <- lapply(files, function(f) {
    d <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (!"iter" %in% names(d) && basename(f) != "lhpars.csv") d$iter <- 1L
    d
  })
  names(out) <- sub("\\.csv$", "", basename(files))
  out
}

readme_text <- function(data, give, answer_key) {
  st <- data$settings
  oem <- st$oem
  c("FLOSS simulated fishery data",
    paste("Created with FLOSS", st$floss_version, "on", format(Sys.Date())),
    "",
    paste("Species:", paste(vapply(st$species, `[[`, "", "species"), collapse = ", ")),
    paste("Fleets:", paste(st$fishery$fleets, collapse = ", ")),
    paste("Years:", paste(range(data$catch$year), collapse = "-")),
    paste("Datasets (iterations):", paste(st$iters, collapse = ", ")),
    "",
    "Files in each iter_xx folder:",
    "  catch.csv   observed catch by fleet (t)",
    "  effort.csv  nominal effort by fleet (fishing days)",
    "  cpue.csv    CPUE by fleet, scaled to mean 1 over the years sampled",
    "  survey.csv  survey biomass index, scaled to mean 1",
    "  lfd.csv     length samples: number of fish per length bin and season (len_mm = lower limit, mm)",
    "  bio.csv     one row per fish: season, len_mm, age (NA = not aged), wt_g, mature (0/1)",
    if (give == "data_lh") "  lhpars.csv  life-history parameters (cm, g, years)",
    "  (optional) catch_aggregate, effort_aggregate, cpue_nominal, lfd_aggregate:",
    "              catch and effort by fleet / all fleets / all species, nominal CPUE",
    "              (kg per fishing day), length samples pooled and raised to the catch",
    "  (optional) floss_flr.rds  FLR objects by species for JABBA3 and FLicc (as_flr())",
    "",
    "Sampling years by fleet:",
    paste0("  ", oem$fleets, ": ", vapply(oem$years, paste, "", collapse = "-")),
    if (answer_key) c("", "answer_key/ holds the true stock dynamics and reference points.",
                      "Remove it before giving the data to trainees."))
}

# Aggregated tables of one iteration (for teaching), without the true values
aggregate_tables <- function(data, i) {
  d <- data
  for (s in c("catch", "effort", "lfd")) d[[s]] <- d[[s]][d[[s]]$iter == i, ]
  lvl <- function(f, levels, ...) do.call(rbind, lapply(levels, function(b) {
    x <- f(d, by = b, ...)
    cbind(level = b, x)
  }))
  drop_true <- function(x) x[, !grepl("_true$|^iter$", names(x)), drop = FALSE]
  list(
    catch_aggregate = drop_true(lvl(catch_table, c("fleet", "aggregate", "combined"))),
    effort_aggregate = drop_true(lvl(effort_table, c("fleet", "aggregate"))),
    cpue_nominal = drop_true(lvl(nominal_cpue, c("fleet", "aggregate", "combined"))),
    lfd_aggregate = drop_true(rbind(
      cbind(raise = "sample", lvl(lfd_aggregate, c("fleet", "aggregate"), raise = "sample")),
      cbind(raise = "catch", lvl(lfd_aggregate, c("fleet", "aggregate"), raise = "catch")))))
}
