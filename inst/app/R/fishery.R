# Step 1: Fishery (Fleets and effort, Effort creep, Effort realisation)

fishery_ui <- function() {
  tabsetPanel(
    id = "fish_tabs",
    tabPanel("Fleets and effort",
      box(width = 12, title = "Fishery",
        fluidRow(
          column(4, selectInput("fish_preset", lab("Fleets", "fleets"),
            choices = stats::setNames(names(fishery_presets),
                                      vapply(fishery_presets, `[[`, "", "label")))),
          column(4, sliderInput("fish_years", "Years (first year = unfished)", min = 1940, max = 2030,
                                value = c(1969, 2025), step = 1, sep = "")),
          column(2, numericInput("fish_seed", "Seed", value = 123, min = 1, step = 1)),
          column(2, numericInput("fish_sigmaE", lab("Effort noise", "sigmaE"), value = 0.15,
                                 min = 0, max = 1, step = 0.05))),
        uiOutput("fish_err"),
        p(class = "note", "The same seed gives the same effort, recruitment and samples every time. All species share one effort history: when a fleet fishes more in a year, it fishes more on every species it catches.")),
      box(width = 12, title = span("Effort trend by fleet", tip("effort_trend")),
        p(class = "note", "Effort starts low (start level) in the first year, rises to the peak by the first change year, moves to the later level by the second change year and stays there. With the one-gear fishery, effort is F/Fmsy: 1 means fishing at Fmsy."),
        uiOutput("fleet_rows"),
        plotOutput("fish_trend", height = 300),
        p(class = "note", "Expected effort without noise. The noise and effort creep are added in the next tabs."))),
    tabPanel("Effort creep",
      box(width = 12, title = span("Effort creep by fleet", tip("creep")),
        p("Fishing efficiency can increase each year (better electronics, gear, skill). Creep raises fishing mortality but is not visible in the nominal effort that fishers report, so CPUE (catch / nominal effort) stays high while the stock declines."),
        uiOutput("creep_rows"),
        plotOutput("fish_creep", height = 340),
        p(class = "note", "Nominal effort (what is reported) and effective effort (nominal x efficiency, what the fish experience)."))),
    tabPanel("Effort realisation",
      box(width = 12, title = "Simulated effort with noise",
        plotOutput("fish_real", height = 380),
        p(class = "note", "Median and 90% interval over iterations, with three single iterations. The number of iterations is set in Operating model.")))
  )
}

fishery_server <- function(input, output, session, rv, r) {

  output$guide_tips <- renderTable(guide_tips_table())

  output$fish_err <- renderUI({
    f <- r$fishery()
    if (inherits(f, "floss_error")) div(class = "warnbox", tags$strong("Check the fishery. "), f$msg)
  })

  # one row per fleet; values start at the preset
  output$fleet_rows <- renderUI({
    p <- preset(input)
    yrs <- input$fish_years
    lapply(seq_along(p$fleets), function(i) {
      d <- p$trends[[i]]
      by <- round(yrs[1] + d$breaks * diff(yrs))
      div(class = "fleet-row", fluidRow(
        column(2, div(class = "fleet-name", p$fleets[i])),
        column(2, numericInput(paste0("fl_e0_", i), "Start level", value = d$e0, min = 0, step = 0.01)),
        column(2, sliderInput(paste0("fl_fhi_", i), "Peak", min = 0, max = 4, value = d$fhi, step = 0.05)),
        column(2, sliderInput(paste0("fl_flo_", i), "Later level", min = 0, max = 4, value = d$flo, step = 0.05)),
        column(4, sliderInput(paste0("fl_brk_", i), "Change years (peak, later level)",
                              min = yrs[1] + 1, max = yrs[2] - 1, value = by, step = 1, sep = ""))))
    })
  })
  outputOptions(output, "fleet_rows", suspendWhenHidden = FALSE)

  output$creep_rows <- renderUI({
    p <- preset(input)
    yrs <- input$fish_years
    lapply(seq_along(p$fleets), function(i) {
      div(class = "fleet-row", fluidRow(
        column(2, div(class = "fleet-name", p$fleets[i])),
        column(5, sliderInput(paste0("fl_creep_", i), "Creep (% per year)", min = 0, max = 5,
                              value = 0, step = 0.25)),
        column(5, sliderInput(paste0("fl_cstart_", i), "Creep starts in", min = yrs[1], max = yrs[2],
                              value = yrs[1], step = 1, sep = ""))))
    })
  })
  outputOptions(output, "creep_rows", suspendWhenHidden = FALSE)

  output$fish_trend <- renderPlot({
    f <- r$fishery()
    validate(need(inherits(f, "floss_fishery"), "Fix the fishery settings above."))
    yrs <- fishery_years(f)
    d <- do.call(rbind, lapply(f$fleets, function(fl)
      data.frame(year = yrs, effort = trend_values(f$trends[[fl]], yrs), fleet = fl)))
    ggplot(d, aes(year, effort, colour = fleet)) + geom_line(linewidth = 1) +
      (if (r$one_gear()) geom_hline(yintercept = 1, linetype = 2, colour = "grey50")) +
      scale_colour_manual(values = app_cols(f$fleets)) +
      labs(x = "Year", y = if (r$one_gear()) "Effort = F / Fmsy" else "Relative effort", colour = NULL) +
      theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  output$fish_creep <- renderPlot({
    f <- r$fishery()
    validate(need(inherits(f, "floss_fishery"), "Fix the fishery settings first."))
    plot_effort(f, nsim = 50, efficiency = TRUE) + theme_bw(base_size = 14)
  })

  output$fish_real <- renderPlot({
    f <- r$fishery()
    validate(need(inherits(f, "floss_fishery"), "Fix the fishery settings first."))
    plot_effort(f, nsim = input$nsim %||% 20) + theme_bw(base_size = 14)
  })
}
