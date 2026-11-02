# Species specification: life history, MSY scale, selectivity and relative
# catchability by fleet.

#' Selectivity-at-length for one fleet
#'
#' Arguments follow [FLRef::sel_la()]: `"logistic"` uses `L50`, `L95`;
#' `"normal"` uses `Lpeak`, `sd`; `"dnormal"` uses `Lpeak`, `sd_left`,
#' `sd_right`. Lengths in cm.
#'
#' @param type `"logistic"`, `"normal"` or `"dnormal"`.
#' @param L50,L95,Lpeak,sd,sd_left,sd_right shape parameters (cm).
#' @return a list of class `floss_sel`.
#' @export
#' @examples
#' selectivity("dnormal", Lpeak = 20, sd_left = 5, sd_right = 15)
selectivity <- function(type = c("logistic", "normal", "dnormal"), L50 = NULL,
                        L95 = NULL, Lpeak = NULL, sd = NULL, sd_left = NULL,
                        sd_right = NULL) {
  type <- match.arg(type)
  need <- switch(type, logistic = c("L50", "L95"), normal = c("Lpeak", "sd"),
                 dnormal = c("Lpeak", "sd_left", "sd_right"))
  pars <- list(L50 = L50, L95 = L95, Lpeak = Lpeak, sd = sd,
               sd_left = sd_left, sd_right = sd_right)
  miss <- need[vapply(pars[need], is.null, TRUE)]
  if (length(miss)) stop(type, " selectivity needs: ", toString(miss))
  structure(c(list(type = type), lapply(pars[need], as.numeric)),
            class = "floss_sel")
}

#' Define a species
#'
#' @param species scientific name.
#' @param common_name common name.
#' @param lh named list or vector of life-history parameters: `linf`, `k`,
#'   `t0` (von Bertalanffy, cm and years), `a`, `b` (W = aL^b, g and cm),
#'   `lm` (length at 50% maturity, cm), `tm` (age at 50% maturity), `tmax`
#'   (maximum age, also the plus group), `M` (natural mortality), `s`
#'   (steepness), `sigR` (recruitment SD) and `rho` (recruitment
#'   autocorrelation, used when `build_om(ar1 = TRUE)`).
#' @param msy target MSY (t) used to scale the stock.
#' @param sel named list of [selectivity()] objects, one per fleet.
#' @param relq named numeric vector of relative catchability by fleet:
#'   apical F of fleet f is `relq[f] x Fmsy x effort[f]`.
#' @param productivity optional label (e.g. "High").
#' @param id file-safe identifier; default from `species`.
#' @return a list of class `floss_species`.
#' @export
species <- function(species, common_name = "", lh, msy, sel, relq,
                    productivity = "", id = NULL) {
  lh <- as.list(lh)
  need <- c("linf", "k", "t0", "a", "b", "lm", "tm", "tmax", "M", "s", "sigR")
  miss <- setdiff(need, names(lh))
  if (length(miss)) stop("missing life-history parameter(s): ", toString(miss))
  if (is.null(lh$rho)) lh$rho <- 0
  lh <- lapply(lh[c(need, "rho")], as.numeric)
  lh$tmax <- as.integer(round(lh$tmax))
  sel <- lapply(sel, function(z) {
    if (!inherits(z, "floss_sel")) z <- do.call(selectivity, z)
    z
  })
  relq <- unlist(relq)
  if (is.null(names(sel)) || is.null(names(relq)))
    stop("'sel' and 'relq' must be named by fleet")
  if (!setequal(names(sel), names(relq)))
    stop("'sel' and 'relq' must name the same fleets")
  out <- structure(list(
    id = if (is.null(id)) species_id(species) else id,
    species = species, common_name = common_name,
    productivity = productivity, version = floss_spec_version(),
    lh = lh, msy = as.numeric(msy),
    sel = sel[names(relq)], relq = as.list(relq)), class = "floss_species")
  validate_species(out)
}

#' @rdname species
#' @param x a `floss_species`.
#' @param fishery optional `floss_fishery`: checks that the species has
#'   selectivity and catchability for every fleet.
#' @export
validate_species <- function(x, fishery = NULL) {
  stopifnot(inherits(x, "floss_species"))
  lh <- x$lh
  if (lh$linf <= 0 || lh$k <= 0) stop("linf and k must be positive")
  if (lh$s <= 0.2 || lh$s > 1) stop("steepness must be in (0.2, 1]")
  if (lh$tmax < 3) stop("tmax must be >= 3")
  if (x$msy <= 0) stop("msy must be positive")
  if (any(unlist(x$relq) < 0)) stop("relq must be >= 0")
  if (!is.null(fishery)) {
    miss <- setdiff(fishery$fleets, names(x$sel))
    if (length(miss)) stop(x$species, ": no selectivity/catchability for ",
                           "fleet(s) ", toString(miss))
  }
  x
}

species_id <- function(x) gsub("[^a-z0-9]+", "_", tolower(trimws(x)))

#' @exportS3Method base::print
print.floss_species <- function(x, ...) {
  cat("FLOSS species: ", x$species,
      if (nzchar(x$common_name)) paste0(" (", x$common_name, ")"), "\n",
      sep = "")
  cat("  MSY =", x$msy, " productivity:", x$productivity, "\n")
  lh <- unlist(x$lh)
  cat("  life history:", paste(names(lh), signif(lh, 3), sep = "=",
                               collapse = ", "), "\n")
  sel <- vapply(x$sel, function(s) paste0(s$type, "(", paste(
    names(s)[-1], unlist(s[-1]), sep = "=", collapse = ", "), ")"), "")
  print(data.frame(fleet = names(x$sel), selectivity = sel,
                   relq = unlist(x$relq), row.names = NULL), row.names = FALSE)
  invisible(x)
}
