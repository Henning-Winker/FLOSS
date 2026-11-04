# The guides ship in inst/guide (always installed); the vignettes include
# them. floss_guide() copies a guide to a folder as a stand-alone R Markdown
# file (and R script), and can render it.

floss_guides <- list(
  `user-guide`   = list(file = "user-guide.Rmd", out = "floss-02-user-guide",
                        title = "FLOSS user guide: four reef species, three fleets, from operating model to saved datasets"),
  teaching       = list(file = "teaching.Rmd", out = "floss-03-teaching",
                        title = "FLOSS teaching scenarios: predict the response, spot the warning flags"),
  `cpue-aggregation` = list(file = "cpue-aggregation.Rmd", out = "floss-04-cpue-aggregation",
                        title = "Training example: summed effort and nominal CPUE hide a shark decline"),
  methods        = list(file = "methods.Rmd", out = "floss-05-methods",
                        title = "Appendix: methods and equations"),
  `how-it-works` = list(file = "how-it-works.md", out = "floss-00-how-it-works",
                        title = "How FLOSS works"),
  `reef-om`      = list(file = "reef-om.Rmd", out = "floss-01-reef-om",
                        title = "FLOSS 1: a multi-fleet reef operating model from saved specs"))

#' Copy, and optionally render, a FLOSS guide
#'
#' The guides are installed with the package (folder `guide`), whether or
#' not the vignettes were built. `floss_guide()` writes a guide to `dir` as a
#' stand-alone R Markdown file that you can edit and knit, plus an R script
#' with its code, and with `render = TRUE` renders it to HTML.
#'
#' @param guide `"user-guide"`, `"teaching"`, `"cpue-aggregation"`, `"methods"` (appendix with the
#'   equations), `"how-it-works"`, `"reef-om"`,
#'   or `"all"`.
#' @param dir folder to write to (created).
#' @param render also render to HTML (needs 'rmarkdown' and pandoc; runs all
#'   the code, which takes a few minutes for the user guide and the teaching
#'   scenarios).
#' @param script also write the R script: the code with the guide text as
#'   `#'` comments ([knitr::purl()] with `documentation = 2`). Fully
#'   commented scripts of the user guide and the teaching scenarios are also
#'   installed: `system.file("guide", "FLOSS_user_guide.R", package = "FLOSS")`.
#' @param open open the HTML after rendering.
#' @param overwrite replace existing files.
#' @return the paths of the files written, invisibly.
#' @export
#' @examples
#' \dontrun{
#' floss_guide("user-guide")                 # floss-02-user-guide.Rmd and .R
#' floss_guide("all", dir = "guides", render = TRUE)
#' }
floss_guide <- function(guide = c("user-guide", "teaching", "cpue-aggregation", "methods", "how-it-works",
                                  "reef-om", "all"),
                        dir = ".", render = FALSE, script = TRUE, open = interactive(),
                        overwrite = FALSE) {
  guide <- match.arg(guide)
  if (guide == "all") {
    out <- unlist(lapply(names(floss_guides), floss_guide, dir = dir, render = render,
                         script = script, open = FALSE, overwrite = overwrite))
    if (render && open) utils::browseURL(file.path(dir, "floss-02-user-guide.html"))
    return(invisible(out))
  }
  g <- floss_guides[[guide]]
  src <- system.file("guide", g$file, package = "FLOSS")
  if (!nzchar(src)) stop("guide file not found: reinstall FLOSS")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rmd <- file.path(dir, paste0(g$out, ".Rmd"))
  if (file.exists(rmd) && !overwrite)
    stop(rmd, " exists; use overwrite = TRUE to replace it")
  writeLines(c("---", paste0("title: \"", g$title, "\""),
               "output:", "  html_document:", "    toc: true", "    toc_float: true", "---", "",
               readLines(src, warn = FALSE)), rmd)
  out <- rmd
  if (script) {
    if (!requireNamespace("knitr", quietly = TRUE)) stop("script = TRUE needs 'knitr'")
    r <- sub("\\.Rmd$", ".R", rmd)
    knitr::purl(rmd, output = r, documentation = 2, quiet = TRUE)   # text kept as #' comments
    out <- c(out, r)
  }
  if (render) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("render = TRUE needs 'rmarkdown'")
    html <- rmarkdown::render(rmd, quiet = TRUE, envir = new.env())
    out <- c(out, html)
    if (open) utils::browseURL(html)
  }
  message("Written: ", paste(basename(out), collapse = ", "))
  invisible(out)
}
