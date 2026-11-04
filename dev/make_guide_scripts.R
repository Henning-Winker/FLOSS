# Regenerates the commented R-script versions of the guides in inst/guide from
# their R Markdown sources (also in inst/guide), as in JABBA3:
#   knitr::purl(documentation = 2) keeps the guide text as #' comments and marks
#   every code chunk "## ---- <chunk name> ----"; chunks the guide shows but
#   does not run (eval = FALSE) are commented out with "## ".
# Run from the package root after editing a guide:  source("dev/make_guide_scripts.R")

guides <- list(
  list(src = "inst/guide/user-guide.Rmd", out = "inst/guide/FLOSS_user_guide.R",
       title = "FLOSS user guide as an R script",
       about = c("Four reef species and three fleets, from the fishery to saved datasets: fishery and effort,",
                 "species, operating models, sampling (seasons, length and biological samples), observations,",
                 "naive aggregation for teaching, true status, export (csv, Excel, FLR objects for JABBA3 and",
                 "FLicc), treatments and an effort-creep teaching scenario."),
       needs = "FLOSS and its FLR dependencies; writexl for Excel export; FLicc and JABBA3 only for the commented-out examples."),
  list(src = "inst/guide/teaching.Rmd", out = "inst/guide/FLOSS_teaching.R",
       title = "FLOSS teaching scenarios as an R script",
       about = c("Teaching cards: predict how the stock responds to a treatment, then check the true dynamics and",
                 "the data, and learn which patterns are warning flags rather than real dynamics."),
       needs = "FLOSS and its FLR dependencies."),
  list(src = "inst/guide/methods.Rmd", out = "inst/guide/FLOSS_methods.R",
       title = "FLOSS appendix: methods and equations, as an R script",
       about = c("How FLOSS simulates a fishery and its data: effort, fishing mortality by gear and species,",
                 "biology and MSY calibration, age-structured dynamics, treatments and sampling (two-stage",
                 "length sampling, autocorrelated CPUE, seasons, ageing error), with small R examples."),
       needs = "FLOSS and its FLR dependencies."),
  list(src = "inst/guide/cpue-aggregation.Rmd", out = "inst/guide/FLOSS_cpue_aggregation.R",
       title = "FLOSS training example as an R script: summed effort and nominal CPUE",
       about = c("Blacktip reef shark in the three-fleet reef fishery: why nominal CPUE from effort summed over",
                 "fleets hides the decline, catchability per fishing day, effective effort and a standardised index."),
       needs = "FLOSS and its FLR dependencies."))

for (g in guides) {
  tmp <- tempfile(fileext = ".R")
  knitr::purl(g$src, output = tmp, documentation = 2, quiet = TRUE)
  body <- readLines(tmp)
  # drop the styling chunk and the knitr set-up chunk (not needed in a script)
  body <- body[!grepl("^## \\.(info|flag|predict)", body)]
  i <- grep("^## ----setup", body)
  if (length(i)) {
    j <- i + which(body[(i + 1):length(body)] == "")[1]
    body <- body[-(i:j)]
  }
  # boxes (::: {.info} ... :::) become plain-text labels
  body <- sub("^#' ::: \\{\\.info\\}$", "#' NOTE:", body)
  body <- sub("^#' ::: \\{\\.flag\\}$", "#' WARNING FLAG:", body)
  body <- sub("^#' ::: \\{\\.predict\\}$", "#' PREDICT FIRST:", body)
  body <- body[!grepl("^#' :::$", body)]
  # no more than one empty comment line in a row
  blank <- body == "#' "
  body <- body[!(blank & c(FALSE, utils::head(blank, -1)))]
  rule <- strrep("=", 117)
  header <- c(
    paste("#", rule),
    paste("#", g$title),
    paste("#", strrep("-", 117)),
    paste("#", g$about),
    "#",
    paste0("# Generated from ", g$src, " with knitr::purl(documentation = 2) (dev/make_guide_scripts.R): the guide text"),
    "# is kept as #' comments, section headings as \"#' #\" lines, and every code chunk is marked",
    "# \"## ---- <chunk name> ----\". Chunks the guide shows but does not run (eval = FALSE) are commented out with \"## \".",
    "#",
    paste("# Requirements:", g$needs),
    "# Install: install.packages(c(\"FLCore\", \"ggplotFL\", \"FLBRP\", \"FLasher\", \"FLife\", \"FLRef\"),",
    "#                           repos = \"https://flr.r-universe.dev\")",
    "#          remotes::install_github(\"Henning-Winker/FLOSS\")",
    "# Run section by section (RStudio: Ctrl+Enter, or Ctrl+Alt+T for the current chunk). Objects are reused in later",
    "# sections, so run from the top. The rendered guide: https://henning-winker.github.io/FLOSS/",
    paste("#", rule), "")
  writeLines(c(header, body), g$out)
  message("written ", g$out, " (", length(body), " lines)")
}
