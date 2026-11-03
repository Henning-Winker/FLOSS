# Step 5: Sampling (Overview, Catch, Effort, CPUE, Survey, Length frequencies,
# Biological samples). Settings per fleet, one row each; press Sample data in
# the left panel to draw the data.

sampling_ui <- function() {
  tagList(
    uiOutput("sm_empty"),
    tabsetPanel(
      id = "sm_tabs",
      tabPanel("Overview",
        box(width = 12, title = "Which years each fleet is sampled",
          fluidRow(column(4, numericInput("samp_n", lab("Datasets to sample", "nsim"), value = 5, min = 1, max = 100, step = 1))),
          p(class = "note", "Each dataset comes from one iteration of the operating model (iterations 1, 2, ...). Choose which to show with Iteration at the top. Default sampling: the most recent years only, longer for fleets with a monitoring history."),
          uiOutput("sm_year_rows")),
        box(width = 12, title = "Data available", plotOutput("sm_avail", height = 340)),
        uiOutput("sm_planted")),
      tabPanel("Catch",
        box(width = 12, title = "Reported catch by fleet",
          p(class = "note", "Observation error (CV) and under-reporting by fleet."),
          uiOutput("ct_rows"),
          fluidRow(column(4, checkboxInput("ct_from_on", "Under-reporting starts in a given year", FALSE)),
                   column(4, conditionalPanel("input.ct_from_on", uiOutput("ct_from_ui"))))),
        box(width = 12, title = "Catch",
          radioButtons("ct_by", NULL, inline = TRUE, choices = c("By fleet" = "fleet",
            "All fleets summed" = "aggregate", "All species and fleets summed" = "combined")),
          plotOutput("sm_catch", height = 380),
          p(class = "note", "Points: reported catch. Dashed (Hide truth unticked): true catch."))),
      tabPanel("Effort",
        box(width = 12, title = "Reported effort by fleet",
          uiOutput("ef_rows"),
          p(class = "note", "Effort is reported in fishing days per year: relative effort (Fishery step) x the fleet's fishing days at relative effort 1, with observation error. It does not include effort creep: that is the point of creep.")),
        box(width = 12, title = "Effort",
          radioButtons("ef_by", NULL, inline = TRUE, choices = c("By fleet" = "fleet", "All fleets summed" = "aggregate")),
          plotOutput("sm_effort", height = 360),
          div(class = "flagbox", tags$strong("Teaching: "), "summing fishing days over fleets treats a day of trap fishing and a day of linefishing as the same effort, although they catch different amounts and sizes of fish."))),
      tabPanel("CPUE",
        box(width = 12, title = "Catch per unit effort by fleet",
          uiOutput("cp_rows"),
          fluidRow(
            column(4, checkboxInput("cp_creep", "CPUE includes effort creep", TRUE)),
            column(4, checkboxInput("cp_qstep_on", lab("Catchability step (warning flag)", "q_step"), FALSE))),
          conditionalPanel("input.cp_qstep_on", uiOutput("cp_qstep_ui")),
          p(class = "note", "CPUE = q x vulnerable biomass^beta x error. Errors are lognormal with autocorrelation rho; high rho mimics non-random sampling (a fleet returning to the same grounds), so residuals come in runs.")),
        box(width = 12, title = "Indices",
          radioButtons("cp_view", NULL, inline = TRUE, choices = c(
            "Index by fleet (as sampled)" = "index", "Nominal CPUE by fleet" = "fleet",
            "Nominal: all fleets" = "aggregate", "Nominal: all species and fleets" = "combined")),
          plotOutput("sm_cpue", height = 400),
          p(class = "note", "Index by fleet: the CPUE series of the sampled years and the survey, each relative to its mean. Nominal CPUE: reported catch / reported fishing days for all years (kg per day, shown relative to its mean), by fleet, summed over fleets, or summed over species and fleets. Dashed (Hide truth unticked): true catch / true effort, or the perfect index; grey: true spawning biomass."),
          div(class = "flagbox", tags$strong("Teaching: "), "nominal CPUE summed over fleets changes when the fleet mix changes, even if the stock does not; summed over species, it can hide the decline of one species behind the others."))),
      tabPanel("Survey",
        box(width = 12, title = "Fishery-independent survey",
          checkboxInput("sv_on", "Include a survey", TRUE),
          conditionalPanel("input.sv_on", uiOutput("sv_ui"))),
        box(width = 12, title = "Survey index", plotOutput("sm_survey", height = 340))),
      tabPanel("Length frequencies",
        box(width = 12, title = "Length sampling (two stages)",
          p(class = "note", "Stage 1: clusters (trips or hauls) are sampled from the catch; fish in one cluster are alike. Stage 2: fish are measured from those clusters. Few clusters give noisy, lumpy length distributions even when many fish are measured."),
          uiOutput("lf_rows"),
          fluidRow(column(4, numericInput("lf_bin", "Length bin (mm, 0 = automatic)", value = 0, min = 0, step = 5)))),
        box(width = 12, title = "Length samples",
          fluidRow(column(4, radioButtons("lf_by", "Fleets", inline = TRUE,
                                          choices = c("By fleet" = "fleet", "All fleets pooled" = "aggregate"))),
                   column(5, radioButtons("lf_raise", "Numbers", inline = TRUE,
                                          choices = c("As sampled" = "sample", "Raised to the catch" = "catch")))),
          fluidRow(column(3, uiOutput("lf_src_ui")),
                   column(3, radioButtons("lf_type", "View", inline = TRUE,
                                          choices = c("Bubbles (all years)" = "bubble", "Histograms" = "bars"))),
                   column(3, conditionalPanel("input.lf_type == 'bars'", uiOutput("lf_years_ui"))),
                   column(3, checkboxInput("lf_coh", "Show year classes (growth lines or ages)", TRUE))),
          plotOutput("sm_lfd", height = 560),
          p(class = "note", "Bubble size: share of fish in each length bin. Real modes move along the growth lines from year to year; modes that jump or vanish are sampling noise. Pooled as sampled, each fleet counts by its sample size; raised to the catch, by its catch (numbers = sample x catch / sample weight)."),
          div(class = "flagbox", tags$strong("Teaching: "), "pooling the samples of all fleets as measured lets the most heavily sampled fleet dominate the length distribution, not the fleet that catches most. Compare pooled and raised.")),
        box(width = 12, title = "Length-based indicators",
          plotOutput("sm_ind", height = 460),
          div(class = "flagbox", tags$strong("Warning flags: "),
              "mean length falling while the index rises; Z from mean length moving without any change in effort (a strong year class); fish above Linf."))),
      tabPanel("Seasons",
        box(width = 12, title = span("Seasonal sampling", tip("seasons")),
          p("Length and biological samples are taken in seasons through the year. Fish grow between seasons and older fish die, so samples from late in the year have larger fish. If the sample is spread unevenly over the seasons, or the sampling season changes, the pooled annual length distribution shifts although the stock has not changed."),
          fluidRow(
            column(2, selectInput("se_n", "Seasons", choices = c(1, 2, 3, 4, 6, 12), selected = 4)),
            column(5, radioButtons("se_type", "Design", choices = c(
              "Balanced: equal samples in every season" = "balanced",
              "Fixed imbalance: unequal shares" = "fixed",
              "Shift: shares change from a given year" = "shift",
              "Opportunistic: shares vary at random each year" = "random"))),
            column(5, uiOutput("se_fleets_ui"))),
          conditionalPanel("input.se_type != 'balanced'", uiOutput("se_share_ui")),
          conditionalPanel("input.se_type == 'shift'", uiOutput("se_shift_ui")),
          conditionalPanel("input.se_type == 'random'",
            sliderInput("se_kappa", lab("How regular (higher = closer to the shares above)", "season_kappa"),
                        min = 1, max = 100, value = 5, step = 1)),
          p(class = "note", "Shares are relative: they are rescaled to sum to 1. The survey keeps its own timing (Survey tab). The stock itself is annual with effort spread evenly through the year: seasons only set when fish are measured, which spreads their lengths around the growth curve."),
          tags$hr(),
          tags$strong("Sampling intensity over the years (same fleets)"),
          fluidRow(
            column(4, checkboxInput("se_nm_on", "Sample size changes from a given year", FALSE)),
            column(4, sliderInput("se_pmiss", "Chance that a year is not sampled", min = 0, max = 0.6, value = 0, step = 0.05))),
          conditionalPanel("input.se_nm_on", uiOutput("se_nm_ui")),
          div(class = "flagbox", tags$strong("Warning flag: "),
              "mean length that steps up or down in the year the sampling season changed, or that tracks the share of samples from one season, is a sampling effect, not a change in growth or mortality.")),
        box(width = 12, title = "Samples and mean length by season",
          plotOutput("sm_seasons", height = 520),
          p(class = "note", "Top: fish measured by season. Bottom: mean length by season (points) and pooled for the year (solid: counts summed; dashed: every sampled season weighted equally). Equal weights correct an uneven spread, but cannot recover a season that was not sampled."))),
      tabPanel("Biological samples",
        box(width = 12, title = "Fish measured, weighed, aged and staged",
          uiOutput("bi_rows"),
          fluidRow(
            column(3, selectInput("bi_scheme", lab("Ageing scheme", "age_scheme"),
                                  choices = c("Random share" = "random", "Length-stratified" = "stratified"))),
            column(3, conditionalPanel("input.bi_scheme == 'random'",
              sliderInput("bi_p_aged", "Share of fish aged", min = 0, max = 1, value = 0.5, step = 0.05)),
              conditionalPanel("input.bi_scheme == 'stratified'",
              numericInput("bi_n_per_bin", "Otoliths per length bin", value = 5, min = 1))),
            column(2, numericInput("bi_age_cv", "Ageing CV", value = 0.1, min = 0, step = 0.05)),
            column(2, numericInput("bi_age_bias", lab("Ageing bias (%)", "age_bias"), value = 0, step = 5)),
            column(2, numericInput("bi_wt_cv", "Weighing CV", value = 0.1, min = 0, step = 0.05)))),
        box(width = 12, title = "Observations",
          fluidRow(
            column(4, selectInput("bi_what", "Show", choices = c("Age and length" = "age_length",
              "Length and weight" = "length_weight", "Maturity at length" = "maturity",
              "Age-length key" = "alk", "Length distribution" = "length"))),
            column(4, uiOutput("bi_src_ui")),
            column(4, radioButtons("bi_by", "Group by", inline = TRUE,
                                   choices = c("Fleet" = "fleet", "All fleets together" = "aggregate", "Season" = "season")))),
          plotOutput("sm_bio", height = 420),
          p(class = "note", "Observations only: nothing is fitted here.")),
        box(width = 12, title = "Catch at age from the length samples and the age-length key",
          plotOutput("sm_caa", height = 380),
          p(class = "note", "Length samples converted to ages with the pooled age-length key of the aged fish. Year classes should appear as diagonals.")))
    )
  )
}

sampling_server <- function(input, output, session, rv, r) {

  spname <- reactive(r$species_used()[[r$shown()]]$species)
  dat <- reactive({
    validate(need(!is.null(rv$data), "Press Sample data in the left panel."))
    validate(need(spname() %in% rv$data$catch$species,
                  "This species has no data yet: build the operating model and sample again."))
    rv$data
  })

  output$sm_empty <- renderUI({
    if (is.null(rv$oms)) empty_msg("Build the operating model first. Then set the sampling here and press Sample data.")
    else if (is.null(rv$data)) empty_msg("Set the sampling below, then press Sample data in the left panel.")
    else if (!r$data_current()) div(class = "warnbox", "Settings changed since the data were sampled. Press Sample data to update.")
  })

  # per-fleet rows: fields is a list of list(id, label, value(i), min, max, step, slider?)
  fleet_rows <- function(fields) {
    p <- preset(input, rv)
    lapply(seq_along(p$fleets), function(i) div(class = "fleet-row", fluidRow(
      column(2, div(class = "fleet-name", p$fleets[i])),
      lapply(fields, function(f) column(f$w %||% 3,
        if (isTRUE(f$slider)) sliderInput(paste0(f$id, "_", fid(p$fleets[i])), f$label, min = f$min, max = f$max,
                                          value = keepv(input, paste0(f$id, "_", fid(p$fleets[i])), f$value(i, p)), step = f$step)
        else numericInput(paste0(f$id, "_", fid(p$fleets[i])), f$label,
                          value = keepv(input, paste0(f$id, "_", fid(p$fleets[i])), f$value(i, p)), min = f$min,
                          max = f$max, step = f$step))))))
  }
  keep <- function(id) outputOptions(output, id, suspendWhenHidden = FALSE)

  output$sm_year_rows <- renderUI({
    req(r$fishery_ok())
    yrs <- r$fishery()$years
    p <- preset(input, rv)
    lapply(seq_along(p$fleets), function(i) div(class = "fleet-row", fluidRow(
      column(2, div(class = "fleet-name", p$fleets[i])),
      column(8, sliderInput(paste0("sm_years_", fid(p$fleets[i])), "Years sampled (catch is reported for all years)",
                            min = yrs[1] + 2, max = yrs[2],
                            value = keepv(input, paste0("sm_years_", fid(p$fleets[i])), c(yrs[2] - p$tail[[i]] + 1, yrs[2])),
                            step = 1, sep = "")))))
  })
  keep("sm_year_rows")

  output$ct_rows <- renderUI(fleet_rows(list(
    list(id = "ct_cv", label = lab("CV", "catch_cv"), value = function(i, p) 0.1, min = 0, max = 1, step = 0.05),
    list(id = "ct_bias", label = lab("Under-reported (%)", "catch_bias"), value = function(i, p) 0,
         min = 0, max = 90, step = 5, slider = TRUE, w = 5))))
  keep("ct_rows")
  output$ct_from_ui <- renderUI({
    y <- r$fishery()$years
    sliderInput("ct_from", "From year", min = y[1], max = y[2], value = y[2] - 10, step = 1, sep = "")
  })

  output$ef_rows <- renderUI(fleet_rows(list(
    list(id = "ef_cv", label = "CV", value = function(i, p) 0.1, min = 0, max = 1, step = 0.05),
    list(id = "ef_days", label = "Fishing days at relative effort 1", value = function(i, p) p$days[[i]],
         min = 1, max = 1e6, step = 100))))
  keep("ef_rows")

  output$cp_rows <- renderUI(fleet_rows(list(
    list(id = "cp_cv", label = "CV", value = function(i, p) 0.2, min = 0, max = 1, step = 0.05, w = 2),
    list(id = "cp_rho", label = lab("Autocorrelation rho", "cpue_rho"), value = function(i, p) p$rho[[i]],
         min = 0, max = 0.95, step = 0.05, slider = TRUE, w = 4),
    list(id = "cp_beta", label = lab("Proportionality beta", "cpue_beta"), value = function(i, p) 1,
         min = 0.2, max = 1.5, step = 0.05, slider = TRUE, w = 4))))
  keep("cp_rows")
  output$cp_qstep_ui <- renderUI({
    y <- r$fishery()$years
    fluidRow(
      column(3, selectInput("cp_qstep_fleet", "Fleet", choices = r$fleets())),
      column(5, sliderInput("cp_qstep_year", "From year", min = y[1], max = y[2], value = y[2] - 6, step = 1, sep = "")),
      column(4, sliderInput("cp_qstep_mult", "Catchability multiplier", min = 0.3, max = 3, value = 1.5, step = 0.05)))
  })

  output$sv_ui <- renderUI({
    y <- r$fishery()$years
    tagList(
      fluidRow(
        column(6, sliderInput("sv_years", "Survey years", min = y[1] + 2, max = y[2], value = c(y[2] - 24, y[2]), step = 1, sep = "")),
        column(3, numericInput("sv_cv", "CV", value = 0.25, min = 0, step = 0.05)),
        column(3, sliderInput("sv_timing", "Timing (fraction of year)", min = 0, max = 0.9, value = 0, step = 0.1))),
      fluidRow(
        column(3, sliderInput("sv_l50", lab("L50 (share of Linf)", "survey_sel"), min = 0.02, max = 0.9, value = 0.2, step = 0.01)),
        column(3, sliderInput("sv_l95", "L95 (share of Linf)", min = 0.05, max = 1, value = 0.35, step = 0.01)),
        column(3, numericInput("sv_ess_age", lab("Length samples: clusters", "ess_age"), value = 50, min = 1)),
        column(3, numericInput("sv_ess_len", lab("Fish measured", "ess_len"), value = 500, min = 1))),
      fluidRow(
        column(3, sliderInput("sv_rho", "Autocorrelation rho", min = 0, max = 0.95, value = 0, step = 0.05)),
        column(3, br(), checkboxInput("sv_qstep_on", lab("Catchability step", "q_step"), FALSE)),
        column(6, conditionalPanel("input.sv_qstep_on", fluidRow(
          column(6, numericInput("sv_qstep_year", "From year", value = y[2] - 8, min = y[1], max = y[2])),
          column(6, numericInput("sv_qstep_mult", "Multiplier", value = 0.6, min = 0.1, step = 0.1)))))),
      p(class = "note", "Survey selectivity is set relative to Linf, so one survey gear works for species of different size. A low L50 gives a recruit survey that sees young fish."))
  })
  outputOptions(output, "sv_ui", suspendWhenHidden = FALSE)

  output$lf_rows <- renderUI(fleet_rows(list(
    list(id = "lf_ess_age", label = lab("Clusters per year", "ess_age"), value = function(i, p) p$ess_age[[i]],
         min = 1, max = 500, step = 1),
    list(id = "lf_ess_len", label = lab("Fish measured per year", "ess_len"), value = function(i, p) p$ess_len[[i]],
         min = 1, max = 10000, step = 50))))
  keep("lf_rows")

  # ---- seasons ----
  output$se_fleets_ui <- renderUI({
    fl <- r$fleets()
    checkboxGroupInput("se_fleets", "Fleets with this design (the others are balanced)",
                       choices = fl, selected = fl)
  })
  keep("se_fleets_ui")
  share_sliders <- function(prefix, n, def) {
    fluidRow(lapply(seq_len(n), function(k) column(max(1, floor(12 / n)),
      sliderInput(paste0(prefix, k), paste0("S", k), min = 0, max = 1, value = def[k], step = 0.05))))
  }
  output$se_share_ui <- renderUI({
    n <- as.integer(input$se_n %||% 4)
    tagList(tags$strong(if (identical(input$se_type, "shift")) "Shares by season, before the shift" else "Shares by season"),
            share_sliders("se_share_", n, c(0.6, 0.2, 0.1, 0.1, rep(0.1, 8))[seq_len(n)]))
  })
  keep("se_share_ui")
  output$se_shift_ui <- renderUI({
    n <- as.integer(input$se_n %||% 4)
    y <- r$fishery()$years
    tagList(fluidRow(
              column(8, sliderInput("se_shift_year", "Shift from year", min = y[1] + 2, max = y[2], value = y[2] - 7, step = 1, sep = "")),
              column(4, sliderInput("se_shift_over", "Gradually, over (years; 0 = at once)", min = 0, max = 10, value = 0, step = 1))),
            tags$strong("Shares by season, from that year"),
            share_sliders("se_share2_", n, rev(c(0.6, 0.2, 0.1, 0.1, rep(0.1, 8))[seq_len(n)])))
  })
  keep("se_shift_ui")
  output$se_nm_ui <- renderUI({
    y <- r$fishery()$years
    fluidRow(
      column(5, sliderInput("se_nm_year", "From year", min = y[1] + 2, max = y[2], value = y[2] - 5, step = 1, sep = "")),
      column(4, sliderInput("se_nm_mult", "Sample size multiplier", min = 0.05, max = 3, value = 0.3, step = 0.05)),
      column(3, sliderInput("se_nm_over", "Over (years)", min = 0, max = 10, value = 0, step = 1)))
  })
  keep("se_nm_ui")

  output$sm_seasons <- renderPlot({
    d <- dat()
    validate(need(!is.null(d$lfd$season) && length(unique(d$lfd$season)) > 1,
                  "These data have one season per year. Choose more seasons and press Sample data."))
    plot_seasons(d, spname(), r$iters()[1]) + theme(text = element_text(size = 14))
  })

  output$bi_rows <- renderUI(fleet_rows(list(
    list(id = "bi_n", label = lab("Fish per year", "bio_n"), value = function(i, p) p$bio_n[[i]],
         min = 0, max = 2000, step = 10))))
  keep("bi_rows")

  # ---- plots ----
  output$sm_avail <- renderPlot(plot_availability(dat(), spname(), r$iters()[1]) + theme(text = element_text(size = 14)))

  output$sm_planted <- renderUI({
    req(r$show_truth())
    oem <- r$oem()
    req(!inherits(oem, "floss_error"))
    a <- planted_artefacts(oem, r$fishery())
    box(width = 12, title = "Planted warning flags (answer key; visible because Hide truth is unticked)",
        if (nrow(a)) tags$ul(lapply(seq_len(nrow(a)), function(i)
          tags$li(tags$strong(a$artefact[i]), " (", a$where[i], "): ", a$detail[i])))
        else p("None: the data are a fair picture of the stock, apart from sampling noise."))
  })

  output$sm_catch <- renderPlot({
    by <- input$ct_by %||% "fleet"
    p <- if (by == "fleet") plot_obs(dat(), spname(), r$iters(), "catch", truth = r$show_truth())
         else plot_aggregate(dat(), spname(), r$iters(), "catch", by, truth = r$show_truth())
    p + theme(text = element_text(size = 14))
  })

  output$sm_effort <- renderPlot({
    if (identical(input$ef_by, "aggregate"))
      return(plot_aggregate(dat(), spname(), r$iters(), "effort", "aggregate", truth = r$show_truth()) +
               theme(text = element_text(size = 14)))
    d <- dat()$effort
    d <- d[d$iter %in% r$iters(), ]
    d$fleet <- factor(d$fleet, levels = r$fleets())
    d$iter <- factor(d$iter)
    multi <- nlevels(d$iter) > 1
    # fleets are in separate panels: with several iterations, colour = iteration
    p <- ggplot(d, aes(year, effort, colour = if (multi) iter else fleet, group = iter)) +
      geom_line(linewidth = 0.5) + geom_point(size = 1.2) +
      facet_wrap(~fleet, scales = "free_y") +
      scale_colour_manual(values = if (multi) app_iter_cols(levels(d$iter)) else app_cols(r$fleets())) +
      labs(x = "Year", y = "Reported effort (fishing days)", colour = "Iteration") + theme_bw(base_size = 14) +
      theme(legend.position = if (multi) "bottom" else "none")
    if (r$show_truth()) p <- p + geom_line(aes(y = effort_true), linetype = 2, linewidth = 0.4)
    p
  })

  output$sm_cpue <- renderPlot({
    v <- input$cp_view %||% "index"
    p <- if (v == "index") plot_obs(dat(), spname(), r$iters(), "index", truth = r$show_truth())
         else plot_aggregate(dat(), spname(), r$iters(), "cpue", v, truth = r$show_truth())
    p + theme(text = element_text(size = 14))
  })

  output$sm_survey <- renderPlot({
    d <- dat()$survey
    validate(need(!is.null(d) && nrow(d) > 0, "No survey in these data."))
    d <- d[d$species == spname() & d$iter %in% r$iters(), ]
    d$iter <- factor(d$iter)
    multi <- nlevels(d$iter) > 1
    cols <- if (multi) app_iter_cols(levels(d$iter)) else stats::setNames(app_cols(r$fleets())[["survey"]], levels(d$iter))
    p <- ggplot(d, aes(year, index, group = iter, colour = iter)) + geom_line(linewidth = 0.5) + geom_point() +
      scale_colour_manual(values = cols) +
      labs(x = "Year", y = "Survey index (relative to mean)", title = spname(), colour = "Iteration") +
      theme_bw(base_size = 14) + theme(legend.position = if (multi) "bottom" else "none")
    if (r$show_truth()) p <- p + geom_line(aes(y = index_true), linetype = 2, linewidth = 0.4)
    p
  })

  output$lf_src_ui <- renderUI({
    s <- unique(dat()$lfd$source)
    selectInput("lf_src", "Samples from", choices = c("All" = "__all__", stats::setNames(s, s)))
  })
  output$lf_years_ui <- renderUI({
    y <- range(dat()$lfd$year)
    sliderInput("lf_years", "Years", min = y[1], max = y[2], value = c(max(y[1], y[2] - 5), y[2]),
                step = 1, sep = "")
  })
  output$sm_lfd <- renderPlot({
    src <- if (identical(input$lf_src %||% "__all__", "__all__")) NULL else input$lf_src
    ty <- input$lf_type %||% "bubble"
    yrs <- if (ty == "bars" && length(input$lf_years) == 2) seq(input$lf_years[1], input$lf_years[2])
    by <- input$lf_by %||% "fleet"; ra <- input$lf_raise %||% "sample"
    lf <- if (by == "fleet" && ra == "sample") NULL else lfd_aggregate(dat(), spname(), by, ra)
    if (!is.null(lf)) src <- NULL
    FLOSS::plot_lfd(dat(), spname(), r$iters()[1], source = src, type = ty, years = yrs,
             cohorts = isTRUE(input$lf_coh), lfd = lf) +
      theme(text = element_text(size = 14))
  })
  output$sm_ind <- renderPlot(plot_indicators(dat(), spname(), r$iters()[1], truth = r$show_truth()))

  output$bi_src_ui <- renderUI({
    s <- unique(dat()$bio$source)
    selectInput("bi_src", "Samples from", choices = c("All" = "__all__", stats::setNames(s, s)))
  })
  bsrc <- reactive(if (identical(input$bi_src %||% "__all__", "__all__")) NULL else input$bi_src)
  output$sm_bio <- renderPlot({
    plot_bio(dat(), spname(), r$iters()[1], what = input$bi_what %||% "age_length", source = bsrc(),
             by = input$bi_by %||% "fleet") +
      theme(text = element_text(size = 14))
  })
  output$sm_caa <- renderPlot({
    src <- bsrc() %||% unique(dat()$lfd$source)[1]
    caa <- tryCatch(caa_obs(dat(), spname(), source = src, iter = r$iters()[1]), error = function(e) e)
    validate(need(!inherits(caa, "error"), paste("Catch at age not available:", if (inherits(caa, "error")) conditionMessage(caa))))
    plot_caa(caa) + labs(subtitle = paste("Length samples from", src)) + theme(text = element_text(size = 14))
  })
}
