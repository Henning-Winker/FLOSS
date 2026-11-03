#' Start the FLOSS app
#'
#' Opens the FLOSS Shiny app: build a fishery and species, the operating
#' model, treatments, sampling and export, step by step in the left panel.
#' The app uses the same functions as the package, so everything done in the
#' app can be repeated in a script.
#'
#' @param ... passed to [shiny::runApp()], e.g. `launch.browser = TRUE` or
#'   `port = 3838`.
#' @return called for its side effect (runs the app).
#' @export
#' @examples
#' if (interactive()) run_app()
run_app <- function(...) {
  need <- c("shiny", "shinydashboard", "commonmark")
  miss <- need[!vapply(need, requireNamespace, TRUE, quietly = TRUE)]
  if (length(miss))
    stop("The app needs: ", toString(miss), ". Install with install.packages(c(",
         paste0("\"", miss, "\"", collapse = ", "), "))")
  dir <- system.file("app", package = "FLOSS")
  if (!nzchar(dir)) stop("App folder not found: reinstall FLOSS")
  shiny::runApp(dir, ...)
}
