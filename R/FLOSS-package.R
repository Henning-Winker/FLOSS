#' FLOSS: FL Operating model and Sampling Simulator
#'
#' Multi-fleet, multi-species operating models built from life history,
#' fleet selectivity, relative catchability and a target MSY, driven by
#' shared fleet effort with optional effort creep, plus simulated fishery
#' data for training.
#'
#' Workflow: [fishery()] -> [species()] (or [load_species()]) ->
#' [build_om()] / [build_oms()] -> sampling (v0.2) -> export (v0.2).
#'
#' @keywords internal
#' @import FLCore FLBRP FLasher
#' @importFrom methods as new
"_PACKAGE"
