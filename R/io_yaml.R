# Read and write specs as YAML, and manage the species library.

#' Write and read FLOSS specs as YAML
#'
#' Species and fishery specs are stored as plain YAML so that they can be
#' edited by hand, shared between trainers and kept under version control.
#'
#' @param x a `floss_species` or `floss_fishery`.
#' @param file path to a `.yaml` file.
#' @return `write_spec()` returns `file` invisibly; `read_spec()` returns the
#'   spec, re-validated through its constructor.
#' @export
#' @examples
#' f <- tempfile(fileext = ".yaml")
#' write_spec(reef_fishery(), f)
#' read_spec(f)
write_spec <- function(x, file) {
  kind <- if (inherits(x, "floss_species")) "species" else
    if (inherits(x, "floss_fishery")) "fishery" else
      stop("x must be a floss_species or floss_fishery")
  lst <- c(list(floss = kind), unclass_deep(x))
  yaml::write_yaml(lst, file, precision = 10)
  invisible(file)
}

#' @rdname write_spec
#' @export
read_spec <- function(file) {
  y <- yaml::read_yaml(file)
  switch(y$floss %||% "",
    species = species(
      species = y$species, common_name = y$common_name %||% "",
      lh = y$lh, msy = y$msy, sel = y$sel, relq = unlist(y$relq),
      productivity = y$productivity %||% "", id = y$id),
    fishery = fishery(
      fleets = unlist(y$fleets), years = seq(y$years[[1]], y$years[[2]]),
      trends = y$trends, sigmaE = unlist(y$sigmaE), creep = unlist(y$creep),
      creep_start = unlist(y$creep_start), name = y$name %||% "fishery",
      seed = y$seed %||% 123L),
    stop(file, " is not a FLOSS spec (missing 'floss: species|fishery')"))
}

#' Species library
#'
#' The library combines the read-only example species shipped with the
#' package and the species saved by the user with [save_species()], which are
#' kept in `tools::R_user_dir("FLOSS", "data")` (override with
#' `options(FLOSS.user_dir = ...)`). A user species with the same id as an
#' example takes precedence.
#'
#' @return `species_library()`: a data frame with one row per species
#'   (`id`, `species`, `common_name`, `msy`, `source`, `file`).
#' @export
#' @examples
#' species_library()
#' load_species("siganus_sutor")
species_library <- function() {
  files <- c(example = list_yaml(example_dir("species")),
             user = list_yaml(user_dir("species")))
  if (!length(files))
    return(data.frame(id = character(), species = character(),
                      common_name = character(), msy = numeric(),
                      source = character(), file = character()))
  src <- sub("[0-9]+$", "", names(files))
  rows <- lapply(seq_along(files), function(i) {
    sp <- read_spec(files[[i]])
    data.frame(id = sp$id, species = sp$species, common_name = sp$common_name,
               msy = sp$msy, source = src[i], file = files[[i]])
  })
  lib <- do.call(rbind, rows)
  # user versions override examples with the same id
  lib <- lib[order(lib$source == "example"), ]
  lib <- lib[!duplicated(lib$id), ]
  lib[order(lib$source, lib$id), , drop = FALSE]
}

#' @rdname species_library
#' @param id species id (see `species_library()$id`).
#' @export
load_species <- function(id) {
  lib <- species_library()
  if (!id %in% lib$id)
    stop("species '", id, "' not in library; available: ", toString(lib$id))
  read_spec(lib$file[lib$id == id])
}

#' @rdname species_library
#' @param x a `floss_species`.
#' @param overwrite overwrite an existing user species with the same id.
#' @export
save_species <- function(x, overwrite = FALSE) {
  stopifnot(inherits(x, "floss_species"))
  dir <- user_dir("species")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file <- file.path(dir, paste0(x$id, ".yaml"))
  if (file.exists(file) && !overwrite)
    stop("user species '", x$id, "' exists; use overwrite = TRUE")
  write_spec(x, file)
}

#' @rdname species_library
#' @export
delete_species <- function(id) {
  file <- file.path(user_dir("species"), paste0(id, ".yaml"))
  if (!file.exists(file)) stop("no user species '", id, "' (examples are read-only)")
  invisible(file.remove(file))
}

#' Example fishery shipped with the package
#' @param name fishery name (file in `inst/extdata/fishery`).
#' @export
load_fishery <- function(name = "reef_3fleets") {
  file <- file.path(example_dir("fishery"), paste0(name, ".yaml"))
  if (!file.exists(file)) stop("no example fishery '", name, "'")
  read_spec(file)
}

# helpers ---------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

unclass_deep <- function(x) {
  if (is.list(x)) lapply(unclass(x), unclass_deep) else x
}

list_yaml <- function(dir) {
  if (!dir.exists(dir)) return(character())
  sort(list.files(dir, pattern = "\\.ya?ml$", full.names = TRUE))
}

example_dir <- function(what) system.file("extdata", what, package = "FLOSS")

user_dir <- function(what) {
  base <- getOption("FLOSS.user_dir", tools::R_user_dir("FLOSS", "data"))
  file.path(base, what)
}
