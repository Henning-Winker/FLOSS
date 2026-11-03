# Step 4: Treatments (Build, Predict, Reveal). A treatment changes one thing in
# the operating model from a given year and is run against an untreated
# baseline with the same random numbers, so every difference comes from it.

tr_types <- c("Fishing effort" = "effort", "Recruitment" = "recruitment",
              "Selectivity (gear change)" = "selectivity",
              "Natural mortality" = "M", "Effort creep" = "creep")

pred_items <- c(ssb = "Spawning biomass", fbar = "Fishing mortality", catch = "Catch",
                rec = "Recruitment", lmean = "Mean length in the catch")

treat_ui <- function() {
  tagList(
    uiOutput("tr_empty"),
    tabsetPanel(
      id = "tr_tabs",
      tabPanel("Build",
        box(width = 12, title = span("Treatment", tip("treatment")),
          p("Change one thing from a given year. The treatment is run for the species shown at the top, against a baseline with exactly the same random numbers (effort noise, recruitment, sampling), so every difference you see comes from the treatment."),
          fluidRow(
            column(4, radioButtons("tr_type", "What changes", choices = tr_types)),
            column(8, uiOutput("tr_pars"))),
          fluidRow(
            column(4, textInput("tr_name", "Name", placeholder = "e.g. halve linefish effort")),
            column(4, sliderInput("tr_nsim", lab("Iterations", "nsim"), min = 5, max = 50, value = 20, step = 5))),
          actionButton("tr_run", "Run treatment", class = "btn-primary"),
          p(class = "note", "Runs the baseline and the treatment and samples data for the first five iterations with the current Sampling settings. Takes a few seconds."),
          uiOutput("tr_done"))),
      tabPanel("Predict",
        box(width = 12, title = "Your prediction",
          uiOutput("tr_summary"),
          div(class = "predictbox", tags$strong("Before you look: "),
            "for each quantity, will it go up, down or stay the same compared with the baseline by the end of the time series, and how soon will the change show?"),
          uiOutput("pred_rows"),
          actionButton("pred_lock", "Lock in my prediction", class = "btn-primary"),
          uiOutput("pred_locked"))),
      tabPanel("Reveal",
        uiOutput("rev_gate"),
        uiOutput("rev_out"))
    )
  )
}

treat_server <- function(input, output, session, rv, r) {

  output$tr_empty <- renderUI({
    if (!r$fishery_ok() || !length(tryCatch(r$species_used(), error = function(e) list())))
      empty_msg("Set the fishery and choose species first.")
  })

  # parameters for the chosen kind of treatment
  output$tr_pars <- renderUI({
    req(r$fishery_ok())
    f <- r$fishery()
    fl <- f$fleets
    y <- f$years
    def <- y[2] - 12
    yr <- function(id, label = "From year", value = def)
      sliderInput(id, label, min = y[1] + 2, max = y[2], value = value, step = 1, sep = "")
    sp <- tryCatch(r$species_used()[[r$shown()]], error = function(e) NULL)
    linf <- sp$lh$linf %||% 60
    switch(input$tr_type %||% "effort",
      effort = tagList(
        selectInput("tr_eff_fleet", "Fleet", choices = fl, selected = fl[length(fl)]),
        yr("tr_eff_from"),
        sliderInput("tr_eff_mult", "Effort multiplier (0 = the fleet stops)", min = 0, max = 3, value = 0.5, step = 0.05)),
      recruitment = tagList(
        radioButtons("tr_rec_type", NULL, inline = TRUE,
          choices = c("One strong (or weak) year class" = "pulse", "Steady trend" = "trend", "Regime shift" = "regime")),
        yr("tr_rec_year", "Year (pulse) or from year (trend, regime)"),
        conditionalPanel("input.tr_rec_type != 'trend'",
          sliderInput("tr_rec_mult", "Recruitment multiplier", min = 0, max = 10, value = 8, step = 0.1)),
        conditionalPanel("input.tr_rec_type == 'trend'",
          sliderInput("tr_rec_rate", "Change per year (%)", min = -20, max = 20, value = -10, step = 1))),
      selectivity = tagList(
        selectInput("tr_sel_fleet", "Fleet", choices = fl, selected = fl[length(fl)]),
        yr("tr_sel_from"),
        p(class = "note", "New S-shaped selectivity, e.g. a minimum size or a larger hook or mesh."),
        fluidRow(
          column(6, sliderInput("tr_sel_l50", "L50 (cm)", min = 1, max = ceiling(linf), value = round(0.6 * linf), step = 0.5)),
          column(6, sliderInput("tr_sel_l95", "L95 (cm)", min = 1, max = ceiling(1.2 * linf), value = round(0.7 * linf), step = 0.5)))),
      M = tagList(
        yr("tr_m_from"),
        sliderInput("tr_m_mult", "Natural mortality multiplier", min = 0.3, max = 3, value = 1.5, step = 0.05),
        p(class = "note", "The reference points stay those of the baseline, so status is judged against the old productivity.")),
      creep = tagList(
        selectInput("tr_cr_fleet", "Fleet", choices = fl, selected = fl[length(fl)]),
        yr("tr_cr_from"),
        sliderInput("tr_cr_rate", "Efficiency increase (% per year)", min = 0, max = 10, value = 3, step = 0.5)))
  })

  # the treatment object from the inputs
  make_tr <- function() {
    ty <- input$tr_type
    nm <- trimws(input$tr_name %||% "")
    tr <- switch(ty,
      effort = treatment(effort = list(fleet = input$tr_eff_fleet, from = input$tr_eff_from, mult = input$tr_eff_mult)),
      recruitment = treatment(recruitment = switch(input$tr_rec_type,
        pulse = list(type = "pulse", years = input$tr_rec_year, mult = input$tr_rec_mult),
        trend = list(type = "trend", from = input$tr_rec_year, rate = input$tr_rec_rate / 100),
        regime = list(type = "regime", from = input$tr_rec_year, mult = input$tr_rec_mult))),
      selectivity = treatment(selectivity = list(fleet = input$tr_sel_fleet, from = input$tr_sel_from,
        sel = selectivity("logistic", L50 = input$tr_sel_l50, L95 = max(input$tr_sel_l95, input$tr_sel_l50 + 0.5)))),
      M = treatment(M = list(from = input$tr_m_from, mult = input$tr_m_mult)),
      creep = treatment(creep = list(fleet = input$tr_cr_fleet, from = input$tr_cr_from, rate = input$tr_cr_rate / 100)))
    tr$name <- if (nzchar(nm)) nm else tr_label(tr)
    tr
  }

  observeEvent(input$tr_run, {
    sps <- r$species_used()
    sp <- sps[[r$shown()]]
    tr <- tryCatch(make_tr(), error = function(e) e)
    if (inherits(tr, "error")) { showNotification(conditionMessage(tr), type = "error"); return() }
    oem <- r$oem()
    if (inherits(oem, "floss_error")) oem <- oem_settings(r$fishery())
    n <- input$tr_nsim %||% 20
    res <- NULL
    withProgress(message = "Running the baseline and the treatment.", value = 0.5, {
      res <- tryCatch(run_treatment(sp, r$fishery(), stats::setNames(list(tr), tr$name), nsim = n,
                                    oem = oem, iters = seq_len(min(5, n)),
                                    seed = input$fish_seed %||% 123, ar1 = isTRUE(input$ar1)),
                      error = function(e) e)
    })
    if (inherits(res, "error")) {
      showNotification(paste("The treatment could not be run:", conditionMessage(res)), type = "error", duration = 15)
      return()
    }
    rv$treat <- list(res = res, tr = tr, sp = sp, from = tr_from(tr))
    rv$prediction <- NULL
    rv$revealed <- FALSE
    updateTabsetPanel(session, "tr_tabs", "Predict")
  })

  output$tr_done <- renderUI({
    req(rv$treat)
    div(class = "infobox", "Treatment run: ", tags$strong(rv$treat$tr$name), " for ", rv$treat$sp$species,
        ". Now make your prediction in the Predict tab.")
  })

  # ---- predict ----
  output$tr_summary <- renderUI({
    if (is.null(rv$treat)) return(empty_msg("Run a treatment first (Build tab)."))
    div(tags$p(tags$strong("Treatment: "), rv$treat$tr$name, tags$br(),
               tags$strong("Species: "), rv$treat$sp$species, tags$br(),
               tags$strong("Starts in: "), rv$treat$from))
  })

  output$pred_rows <- renderUI({
    req(rv$treat)
    lapply(names(pred_items), function(k) fluidRow(
      column(3, div(style = "padding-top:8px; font-weight:600;", pred_items[[k]])),
      column(5, radioButtons(paste0("pred_", k), NULL, inline = TRUE, selected = character(0),
                             choices = c("Up", "Down", "No change"))),
      column(4, selectInput(paste0("pred_when_", k), NULL,
                            choices = c("When?" = "", "Within a year", "After 2 to 4 years", "After 5 or more years", "Never")))))
  })

  observeEvent(input$pred_lock, {
    req(rv$treat)
    p <- lapply(names(pred_items), function(k)
      data.frame(item = k, direction = input[[paste0("pred_", k)]] %||% "", when = input[[paste0("pred_when_", k)]] %||% ""))
    rv$prediction <- do.call(rbind, p)
    rv$revealed <- TRUE
    updateTabsetPanel(session, "tr_tabs", "Reveal")
  })

  output$pred_locked <- renderUI({
    if (!is.null(rv$prediction)) div(class = "infobox", "Prediction locked. See the Reveal tab.")
  })

  # ---- reveal ----
  output$rev_gate <- renderUI({
    if (is.null(rv$treat)) return(empty_msg("Run a treatment first (Build tab)."))
    if (isTRUE(rv$revealed)) return(NULL)
    div(class = "predictbox", "Make and lock your prediction first (Predict tab). ",
        actionLink("rev_skip", "Or show the result without a prediction."))
  })
  observeEvent(input$rev_skip, rv$revealed <- TRUE)

  output$rev_out <- renderUI({
    req(rv$treat, isTRUE(rv$revealed))
    fl <- rv$treat$res$oms$baseline[[1]]$fishery$fleets
    lf <- rv$treat$res$data$baseline$lfd
    src <- unique(lf$source)
    # default: the treated fleet if it has length samples, else the longest series
    tf <- c(vapply(rv$treat$tr$effort, `[[`, "", "fleet"), vapply(rv$treat$tr$selectivity, `[[`, "", "fleet"),
            vapply(rv$treat$tr$creep, `[[`, "", "fleet"))
    ny <- vapply(src, function(s) length(unique(lf$year[lf$source == s])), 0)
    src0 <- if (length(tf) && tf[1] %in% src) tf[1] else src[which.max(ny)]
    tagList(
      box(width = 12, title = "Your prediction and what happened",
          tableOutput("rev_table"),
          p(class = "note", "What happened: median over iterations of the last five years, treatment compared with baseline. A change of less than 5% (2% for mean length) counts as no change. When: the first year after the start with a change of that size.")),
      box(width = 12, title = "True stock: baseline and treatment",
          plotOutput("rev_plot", height = 460),
          p(class = "note", "Dashed black: baseline. Both runs share every random number, so the lines are identical until the treatment starts.")),
      box(width = 12, title = "What the length data show",
          fluidRow(column(4, selectInput("rev_src", "Length samples from", choices = src, selected = src0)),
                   column(4, radioButtons("rev_stat", NULL, inline = TRUE,
                                          choices = c("Median of five datasets" = "median", "Each dataset" = "all")))),
          plotOutput("rev_li", height = 460),
          div(class = "flagbox", tags$strong("Look for: "),
              "does mean length or Z from mean length move in the same direction as the true change? Recruitment changes move length indicators even when fishing has not changed, and a strong year class passing through the samples can look like a change in mortality."),
          p(class = "note", "Z from mean length assumes equilibrium (constant recruitment and mortality) and flat-topped selectivity. The dashed line is the true Z.")))
  })

  output$rev_plot <- renderPlot(plot_treatment(rv$treat$res) + theme_bw(base_size = 14) + theme(legend.position = "bottom"))

  li <- reactive({
    req(rv$treat)
    length_indicators(rv$treat$res, rv$treat$sp$species, iter = rv$treat$res$data$baseline$settings$iters)
  })
  output$rev_li <- renderPlot({
    req(input$rev_src)
    plot_length_indicators(li(), source = input$rev_src, stat = input$rev_stat %||% "median")
  })

  output$rev_table <- renderTable({
    req(rv$treat)
    obs <- tr_response(rv$treat, li(), input$rev_src)
    p <- rv$prediction
    out <- data.frame(Quantity = unname(pred_items[obs$item]),
                      `Your prediction` = if (is.null(p)) "" else paste(p$direction, ifelse(nzchar(p$when), paste0("(", p$when, ")"), "")),
                      `What happened` = obs$direction, `Change at the end` = obs$change,
                      `When` = obs$when, check.names = FALSE)
    if (!is.null(p)) out$Match <- ifelse(!nzchar(p$direction), "", ifelse(p$direction == obs$direction, "yes", "no"))
    out
  })
}

# ---- helpers for treatments --------------------------------------------------

tr_from <- function(tr) {
  min(c(vapply(tr$effort, `[[`, 0, "from"), vapply(tr$selectivity, `[[`, 0, "from"),
        vapply(tr$creep, `[[`, 0, "from"), tr$M$from,
        unlist(lapply(tr$recruitment, function(x) x$from %||% min(x$years)))))
}

tr_label <- function(tr) {
  if (length(tr$effort)) { e <- tr$effort[[1]]; return(paste0(e$fleet, " effort x", e$mult, " from ", e$from)) }
  if (length(tr$recruitment)) { x <- tr$recruitment[[1]]
    return(switch(x$type, pulse = paste0("year class ", x$years, " x", x$mult),
                  trend = paste0("recruitment ", x$rate * 100, "% per year from ", x$from),
                  regime = paste0("recruitment x", x$mult, " from ", x$from))) }
  if (length(tr$selectivity)) { s <- tr$selectivity[[1]]; return(paste0(s$fleet, " L50 ", s$sel$L50, " cm from ", s$from)) }
  if (length(tr$creep)) { c <- tr$creep[[1]]; return(paste0(c$fleet, " creep ", c$rate * 100, "% per year from ", c$from)) }
  if (!is.null(tr$M)) return(paste0("M x", tr$M$mult, " from ", tr$M$from))
  "treatment"
}

# Direction, size and timing of the response (treatment vs baseline)
tr_response <- function(x, li, src = NULL) {
  b <- x$res$oms$baseline[[1]]$stock
  t <- x$res$oms[[2]][[1]]$stock
  yrs <- as.numeric(dimnames(b)$year)
  from <- x$from
  med <- function(q) apply(matrix(c(q), nrow = length(yrs)), 1, stats::median)
  one <- function(qb, qt) {
    ratio <- med(qt) / med(qb)
    tail_r <- mean(utils::tail(ratio, 5))
    after <- yrs >= from
    hit <- which(after & abs(ratio - 1) >= 0.05)
    list(r = tail_r, when = if (length(hit)) yrs[hit[1]] - from else NA)
  }
  res <- list(ssb = one(ssb(b), ssb(t)), fbar = one(fbar(b), fbar(t)),
              catch = one(catch(b), catch(t)), rec = one(rec(b), rec(t)))
  # mean length in the samples of one source (median over datasets)
  l <- li[li$source == (src %||% li$source[1]), ]
  lm <- stats::aggregate(Lmean ~ run + year, l, stats::median)
  lb <- lm[lm$run == "baseline", ]; lt <- lm[lm$run != "baseline", ]
  lt <- lt[match(lb$year, lt$year), ]
  ratio <- lt$Lmean / lb$Lmean
  hit <- which(lb$year >= from & abs(ratio - 1) >= 0.02)
  res$lmean <- list(r = mean(utils::tail(ratio, 5)),
                    when = if (length(hit)) lb$year[hit[1]] - from else NA, thr = 0.02)
  dir <- function(z, thr = 0.05) if (abs(z$r - 1) < thr) "No change" else if (z$r > 1) "Up" else "Down"
  when <- function(z) if (is.na(z$when)) "Never" else if (z$when <= 1) "Within a year" else
    if (z$when <= 4) "After 2 to 4 years" else "After 5 or more years"
  data.frame(item = names(res),
             direction = vapply(res, function(z) dir(z, z$thr %||% 0.05), ""),
             change = sprintf("%+.0f%%", 100 * (vapply(res, `[[`, 0, "r") - 1)),
             when = vapply(res, when, ""))
}
