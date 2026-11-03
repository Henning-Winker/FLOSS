# Step 1: Fishery (Fleets and effort, Effort creep, Effort realisation)

fishery_ui <- function() {
  tabsetPanel(
    id = "fish_tabs",
    tabPanel("Fleets and effort",
      box(width = 12, title = "Fishery",
        fluidRow(
          column(4, selectInput("fish_preset", lab("Fleets", "fleets"),
            choices = c(stats::setNames(names(fishery_presets),
                                        vapply(fishery_presets, `[[`, "", "label")),
                        "New fishery: add your own gears" = "custom"))),
          column(4, sliderInput("fish_years", "Years (first year = unfished)", min = 1940, max = 2030,
                                value = c(1969, 2025), step = 1, sep = "")),
          column(2, numericInput("fish_seed", "Seed", value = 123, min = 1, step = 1)),
          column(2, numericInput("fish_sigmaE", lab("Effort noise", "sigmaE"), value = 0.15,
                                 min = 0, max = 1, step = 0.05))),
        conditionalPanel("input.fish_preset == 'custom'",
          div(class = "infobox", tags$strong("New fishery. "),
              "Start with no gears: type a gear name and press Add gear, then set its effort below. Add as many gears as you need. Each species gets a default selectivity for a new gear (S-shaped at the length at maturity), which you set in Species, Selectivity and Exploitation."),
          fluidRow(
            column(4, textInput("gear_name", "Gear name", placeholder = "e.g. handline, trawl, beach seine")),
            column(2, br(), actionButton("gear_add", "Add gear", class = "btn-primary")),
            column(4, uiOutput("gear_rm_ui")),
            column(2, br(), uiOutput("gear_rm_btn")))),
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

  # ---- new fishery: add and remove gears ----
  observeEvent(input$gear_add, {
    nm <- trimws(input$gear_name %||% "")
    if (!nzchar(nm)) { showNotification("Type a gear name first.", type = "warning"); return() }
    if (tolower(nm) %in% c("survey", "all fleets") || fid(nm) %in% fid(rv$gears)) {
      showNotification(paste0("'", nm, "' is already used. Choose another name."), type = "warning"); return()
    }
    rv$gears <- c(rv$gears, nm)
    updateTextInput(session, "gear_name", value = "")
  })
  output$gear_rm_ui <- renderUI({
    if (!length(rv$gears)) return(NULL)
    selectInput("gear_rm", "Remove a gear", choices = rv$gears)
  })
  output$gear_rm_btn <- renderUI(if (length(rv$gears)) actionButton("gear_rm_go", "Remove"))
  observeEvent(input$gear_rm_go, rv$gears <- setdiff(rv$gears, input$gear_rm))

  output$fish_err <- renderUI({
    f <- r$fishery()
    if (inherits(f, "floss_error")) div(class = "warnbox", tags$strong("Check the fishery. "), f$msg)
  })

  # one row per fleet; values start at the preset
  output$fleet_rows <- renderUI({
    p <- preset(input, rv)
    yrs <- input$fish_years
    if (!length(p$fleets)) return(div(class = "infobox", "No gears yet: add a gear above."))
    lapply(seq_along(p$fleets), function(i) {
      d <- p$trends[[i]]
      k <- fid(p$fleets[i])
      by <- round(yrs[1] + d$breaks * diff(yrs))
      by <- pmin(pmax(keepv(input, paste0("fl_brk_", k), by), yrs[1] + 1), yrs[2] - 1)
      div(class = "fleet-row", fluidRow(
        column(2, div(class = "fleet-name", p$fleets[i])),
        column(2, numericInput(paste0("fl_e0_", k), "Start level", value = keepv(input, paste0("fl_e0_", k), d$e0), min = 0, step = 0.01)),
        column(2, sliderInput(paste0("fl_fhi_", k), "Peak", min = 0, max = 4, value = keepv(input, paste0("fl_fhi_", k), d$fhi), step = 0.05)),
        column(2, sliderInput(paste0("fl_flo_", k), "Later level", min = 0, max = 4, value = keepv(input, paste0("fl_flo_", k), d$flo), step = 0.05)),
        column(4, sliderInput(paste0("fl_brk_", k), "Change years (peak, later level)",
                              min = yrs[1] + 1, max = yrs[2] - 1, value = by, step = 1, sep = ""))))
    })
  })
  outputOptions(output, "fleet_rows", suspendWhenHidden = FALSE)

  output$creep_rows <- renderUI({
    p <- preset(input, rv)
    yrs <- input$fish_years
    if (!length(p$fleets)) return(div(class = "infobox", "No gears yet: add a gear in Fleets and effort."))
    lapply(seq_along(p$fleets), function(i) {
      k <- fid(p$fleets[i])
      div(class = "fleet-row", fluidRow(
        column(2, div(class = "fleet-name", p$fleets[i])),
        column(5, sliderInput(paste0("fl_creep_", k), "Creep (% per year)", min = 0, max = 5,
                              value = keepv(input, paste0("fl_creep_", k), 0), step = 0.25)),
        column(5, sliderInput(paste0("fl_cstart_", k), "Creep starts in", min = yrs[1], max = yrs[2],
                              value = keepv(input, paste0("fl_cstart_", k), yrs[1]), step = 1, sep = ""))))
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
