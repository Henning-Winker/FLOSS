# FLOSS Shiny app: entry point. Start it with FLOSS::run_app().
# Shiny loads every file in the R/ folder next to this file before it runs,
# so the functions defined there (one file per step) are available here.

library(shiny)
library(shinydashboard)
library(ggplot2)
library(FLCore)
library(FLOSS)

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "FLOSS", titleWidth = 280),
  app_sidebar(),
  app_body()
)

server <- function(input, output, session) {
  # rv holds the app's shared state: the species being worked on, the built
  # operating models, the sampled data and the treatment runs
  rv <- reactiveValues(sp = initial_species(), oms = NULL, om_snap = NULL,
                       data = NULL, data_snap = NULL, treat = NULL,
                       prediction = NULL, revealed = FALSE, applied = 0,
                       gears = character(0),
                       # survey selectivity, lengths as shares of Linf (shown in cm)
                       sv = list(L50 = 0.2, L95 = 0.3, Lpeak = 0.3, sd = 0.1, sd_left = 0.08, sd_right = 0.2))

  # r holds the reactive objects (fishery, species, sampling settings, ...)
  # used by several steps
  r <- make_reactives(input, rv)

  run_server(input, output, session, rv, r)
  fishery_server(input, output, session, rv, r)
  species_server(input, output, session, rv, r)
  om_server(input, output, session, rv, r)
  treat_server(input, output, session, rv, r)
  sampling_server(input, output, session, rv, r)
  export_server(input, output, session, rv, r)
  truth_server(input, output, session, rv, r)
}

shinyApp(ui, server)
