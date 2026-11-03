# Step 2: Species (Library, Life history, Selectivity, Exploitation)

lh_fields <- list(
  linf = list("Linf (cm) *", "linf"), k = list("k (per year) *", "k"),
  t0 = list("t0 (years)", NULL), M = list("M (per year) *", "M"),
  a = list("a (weight-length, g/cm) *", NULL), b = list("b (weight-length) *", NULL),
  lm = list("Length at 50% maturity (cm) *", NULL), tm = list("Age at 50% maturity", NULL),
  tmax = list("Maximum age (plus group) *", NULL), s = list("Steepness *", "steepness"),
  sigR = list("Recruitment SD *", "sigR"), rho = list("Recruitment autocorrelation", "ar1"))

species_ui <- function() {
  lib <- species_library()
  tagList(
    fluidRow(
      column(4, selectInput("sp_edit", "Species to view and edit",
                            choices = stats::setNames(lib$id, lib$species))),
      column(8, div(style = "margin-top:25px; display:flex; gap:10px; align-items:center;",
        actionButton("sp_apply", "Apply changes", class = "btn-primary"),
        actionButton("sp_save", "Save to my library"),
        uiOutput("sp_status", inline = TRUE)))),
    tabsetPanel(
      id = "sp_tabs",
      tabPanel("Library",
        box(width = 12, title = "Species in the operating model",
          checkboxGroupInput("sp_use", NULL, inline = TRUE,
            choices = stats::setNames(lib$id, lib$species), selected = lib$id),
          p(class = "note", "Tick the species to simulate. They all share the fishery's effort. Species can be added one at a time: species already built are reused when you build again.")),
        box(width = 12, title = "Species library",
          tableOutput("sp_table"),
          p(class = "note", "The four examples are reef species with life histories from FishLife. Species you save go to your own library, and a saved species replaces an example with the same name.")),
        box(width = 12, title = "Add a species",
          fluidRow(
            column(4, textInput("sp_new_name", "Scientific name", placeholder = "e.g. Epinephelus fuscoguttatus")),
            column(3, textInput("sp_new_common", "Common name", placeholder = "e.g. Brown-marbled grouper")),
            column(3, selectInput("sp_new_from", "Start from", choices = stats::setNames(lib$id, lib$species))),
            column(2, br(), actionButton("sp_new", "Add"))),
          p(class = "note", "The new species starts as a copy. Change its values in the next tabs, press Apply changes, and Save to keep it for later sessions."),
          uiOutput("sp_delete_ui"))),
      tabPanel("Life history",
        box(width = 12, title = "Biological parameters",
          p(class = "note", "Fields marked * are required. Lengths in cm, weight in g."),
          uiOutput("lh_inputs"),
          uiOutput("lh_warn")),
        fluidRow(
          box(width = 4, title = "Growth", plotOutput("sp_growth", height = 260)),
          box(width = 4, title = "Maturity at length", plotOutput("sp_mat", height = 260)),
          box(width = 4, title = "Weight at length", plotOutput("sp_wt", height = 260)))),
      tabPanel("Selectivity",
        box(width = 12, title = span("Selectivity by fleet", tip("selectivity")),
          p(class = "note", "S-shaped: the gear keeps every fish above a size. Dome-shaped: large fish escape (e.g. gillnet mesh, traps, fishing grounds). Lengths in cm."),
          uiOutput("sel_rows"),
          plotOutput("sp_sel", height = 320),
          p(class = "note", "Coloured: selectivity by fleet. Dashed: maturity. Where a fleet's curve sits left of maturity, it catches fish before they spawn.")),
        box(width = 12, title = "Yield and spawning biomass under each fleet's selectivity",
          p("The same stock fished with each fleet's selectivity alone: equilibrium yield and spawning biomass against fishing mortality. Catching fish young lowers the maximum yield (growth overfishing); dome shapes leave a refuge of large spawners."),
          actionButton("sp_eq_go", "Draw"),
          plotOutput("sp_eq", height = 360))),
      tabPanel("Exploitation",
        box(width = 12, title = "Scale and catchability",
          uiOutput("rq_rows"),
          p(class = "note", "Apical F of a fleet = relative catchability x Fmsy x effort. With relative catchability 1 and effort 1, the fleet alone fishes at Fmsy (for its own selectivity, approximately)."),
          uiOutput("rq_onegear")))
    )
  )
}

species_server <- function(input, output, session, rv, r) {

  cur <- reactive({ req(input$sp_edit %in% names(rv$sp)); rv$sp[[input$sp_edit]] })

  # the species being edited, adapted to the fleets of the fishery
  cur_adapted <- reactive({
    req(r$fishery_ok())
    adapt_species(cur(), r$fishery(), r$one_gear())
  })

  # ---- read the inputs back into a species ----
  collect <- function() {
    base <- cur()
    lh <- lapply(stats::setNames(names(lh_fields), names(lh_fields)),
                 function(n) input[[paste0("lh_", n)]] %||% base$lh[[n]])
    fl <- r$fleets()
    sel <- base$sel
    relq <- base$relq
    for (i in seq_along(fl)) {
      ty <- input[[paste0("sel_type_", fid(fl[i]))]]
      if (is.null(ty)) next
      g <- function(n) input[[paste0("sel_", n, "_", fid(fl[i]))]]
      sel[[fl[i]]] <- switch(ty,
        logistic = selectivity("logistic", L50 = g("L50"), L95 = max(g("L95"), g("L50") + 0.5)),
        normal   = selectivity("normal", Lpeak = g("Lpeak"), sd = g("sd")),
        dnormal  = selectivity("dnormal", Lpeak = g("Lpeak_d"), sd_left = g("sdl"), sd_right = g("sdr")))
      relq[[fl[i]]] <- input[[paste0("rq_", fid(fl[i]))]] %||% relq[[fl[i]]] %||% 0.5
    }
    species(base$species, base$common_name, lh = lh, msy = input$lh_msy %||% base$msy,
            sel = sel, relq = unlist(relq[names(sel)]), productivity = base$productivity, id = base$id)
  }
  collected <- reactive(tryCatch(collect(), error = function(e) e))

  output$sp_status <- renderUI({
    x <- collected()
    if (inherits(x, "error")) return(tags$span(style = "color:#7a4b00; font-weight:600;", "Check the values"))
    edited <- !isTRUE(all.equal(unclass(x), unclass(cur()), tolerance = 1e-6))
    tags$span(style = paste0("font-weight:600; color:", if (edited) "#7a4b00" else "#1c6b3d"),
              if (edited) "Edited: press Apply changes" else "Matches the species in use")
  })

  observeEvent(input$sp_apply, {
    x <- collected()
    if (inherits(x, "error")) {
      showNotification(paste("Check the values:", conditionMessage(x)), type = "error", duration = 10)
      return()
    }
    rv$sp[[x$id]] <- x
    rv$applied <- rv$applied + 1
    showNotification(paste(x$species, "updated. Build the operating model again to use it."), type = "message")
  })

  observeEvent(input$sp_save, {
    x <- collected()
    if (inherits(x, "error")) { showNotification(conditionMessage(x), type = "error"); return() }
    rv$sp[[x$id]] <- x
    res <- tryCatch(save_species(x, overwrite = TRUE), error = function(e) e)
    if (inherits(res, "error")) showNotification(conditionMessage(res), type = "error")
    else showNotification(paste(x$species, "saved to your library."), type = "message")
    rv$applied <- rv$applied + 1
  })

  # ---- library ----
  output$sp_table <- renderTable({
    lib <- species_library()
    do.call(rbind, lapply(rv$sp, function(s) data.frame(
      Species = s$species, `Common name` = s$common_name, Productivity = s$productivity,
      `MSY (t)` = s$msy, `Linf (cm)` = s$lh$linf, k = s$lh$k, M = round(s$lh$M, 3),
      Steepness = s$lh$s, Library = if (s$id %in% lib$id[lib$source == "user"]) "yours" else
        if (s$id %in% lib$id) "example" else "this session only",
      check.names = FALSE)))
  })

  update_lists <- function(select = NULL) {
    ch <- stats::setNames(names(rv$sp), vapply(rv$sp, `[[`, "", "species"))
    updateCheckboxGroupInput(session, "sp_use", choices = ch,
                             selected = unique(c(input$sp_use, select)), inline = TRUE)
    updateSelectInput(session, "sp_edit", choices = ch, selected = select %||% input$sp_edit)
    updateSelectInput(session, "sp_new_from", choices = ch)
  }

  observeEvent(input$sp_new, {
    nm <- trimws(input$sp_new_name)
    if (!nzchar(nm)) { showNotification("Enter a scientific name.", type = "warning"); return() }
    from <- rv$sp[[input$sp_new_from]]
    x <- species(nm, trimws(input$sp_new_common), lh = from$lh, msy = from$msy,
                 sel = from$sel, relq = unlist(from$relq), productivity = from$productivity)
    if (x$id %in% names(rv$sp)) { showNotification("That species already exists.", type = "warning"); return() }
    rv$sp[[x$id]] <- x
    update_lists(select = x$id)
    updateTextInput(session, "sp_new_name", value = "")
    updateTextInput(session, "sp_new_common", value = "")
    updateTabsetPanel(session, "sp_tabs", "Life history")
  })

  output$sp_delete_ui <- renderUI({
    rv$applied
    lib <- species_library()
    mine <- lib$id[lib$source == "user"]
    if (!length(mine)) return(NULL)
    fluidRow(
      column(4, selectInput("sp_del", "Remove from your library", choices = stats::setNames(mine, lib$species[lib$source == "user"]))),
      column(2, br(), actionButton("sp_del_go", "Remove")))
  })
  observeEvent(input$sp_del_go, {
    id <- input$sp_del
    delete_species(id)
    lib <- species_library()
    if (id %in% lib$id) rv$sp[[id]] <- load_species(id) else rv$sp[[id]] <- NULL
    rv$applied <- rv$applied + 1
    update_lists()
  })

  # ---- life history ----
  output$lh_inputs <- renderUI({
    rv$applied
    s <- cur()
    boxes <- lapply(names(lh_fields), function(n) {
      f <- lh_fields[[n]]
      column(2, numericInput(paste0("lh_", n), lab(f[[1]], f[[2]]), value = s$lh[[n]],
                             step = if (n %in% c("tmax")) 1 else NA))
    })
    tagList(fluidRow(boxes[1:6]), fluidRow(boxes[7:12]),
            fluidRow(column(3, numericInput("lh_msy", lab("Target MSY (t) *", "msy"), value = s$msy, min = 0))))
  })
  outputOptions(output, "lh_inputs", suspendWhenHidden = FALSE)

  lh_now <- reactive({
    s <- cur()
    v <- lapply(stats::setNames(names(lh_fields), names(lh_fields)),
                function(n) input[[paste0("lh_", n)]] %||% s$lh[[n]])
    v
  })

  output$lh_warn <- renderUI({
    lh <- lh_now()
    w <- character(0)
    if (any(vapply(lh, function(v) is.null(v) || is.na(v), TRUE))) w <- c(w, "Some values are blank.")
    else {
      if (lh$lm >= lh$linf) w <- c(w, "Length at maturity should be below Linf.")
      if (lh$b < 2 || lh$b > 4) w <- c(w, "b is usually between 2.5 and 3.5. Check the value and its units.")
      if (lh$M / lh$k < 0.5 || lh$M / lh$k > 3) w <- c(w, paste0("M/k = ", round(lh$M / lh$k, 2), ", outside the usual 0.5 to 3. Check M and k."))
      if (lh$s <= 0.2 || lh$s > 1) w <- c(w, "Steepness must be above 0.2 and at most 1.")
      if (lh$tmax < 3) w <- c(w, "Maximum age must be 3 or more.")
    }
    if (length(w)) div(class = "warnbox", tags$strong("Check these values. "), tags$ul(lapply(w, tags$li)))
  })

  vb <- function(lh, a) lh$linf * (1 - exp(-lh$k * (a - lh$t0)))
  okl <- reactive({ lh <- lh_now(); req(!any(vapply(lh, function(v) is.null(v) || is.na(v), TRUE))); lh })

  output$sp_growth <- renderPlot({
    lh <- okl()
    a <- seq(0, lh$tmax, length.out = 200)
    ggplot(data.frame(age = a, len = vb(lh, a)), aes(age, len)) + geom_line(linewidth = 1, colour = "#1874cd") +
      geom_hline(yintercept = lh$linf, linetype = 2, colour = "grey50") +
      labs(x = "Age (years)", y = "Length (cm)") + theme_bw(base_size = 13)
  })
  output$sp_mat <- renderPlot({
    lh <- okl()
    L <- seq(0, lh$linf * 1.1, length.out = 200)
    m <- 1 / (1 + exp(-log(19) * (L - lh$lm) / (0.15 * lh$lm)))
    ggplot(data.frame(L, m), aes(L, m)) + geom_line(linewidth = 1, colour = "#1874cd") +
      geom_vline(xintercept = lh$lm, linetype = 2, colour = "grey50") +
      labs(x = "Length (cm)", y = "Proportion mature") + theme_bw(base_size = 13)
  })
  output$sp_wt <- renderPlot({
    lh <- okl()
    L <- seq(0, lh$linf * 1.1, length.out = 200)
    ggplot(data.frame(L, w = lh$a * L^lh$b / 1000), aes(L, w)) + geom_line(linewidth = 1, colour = "#1874cd") +
      labs(x = "Length (cm)", y = "Weight (kg)") + theme_bw(base_size = 13)
  })

  # ---- selectivity: one row per fleet ----
  output$sel_rows <- renderUI({
    rv$applied
    s <- cur_adapted()
    fl <- r$fleets()
    Lmax <- ceiling(s$lh$linf * 1.3)
    lapply(seq_along(fl), function(i) {
      z <- s$sel[[fl[i]]]
      v <- function(n, def) round(z[[n]] %||% def, 1)
      lm <- s$lh$lm
      id <- function(n) paste0("sel_", n, "_", fid(fl[i]))
      sl <- function(n, label, def) sliderInput(id(n), label, min = 0, max = Lmax, value = v(n, def), step = 0.5)
      cond <- function(t) paste0("input.sel_type_", fid(fl[i]), " == '", t, "'")
      div(class = "fleet-row", fluidRow(
        column(2, div(class = "fleet-name", fl[i])),
        column(3, selectInput(paste0("sel_type_", fid(fl[i])), "Shape", choices = sel_types, selected = z$type)),
        column(7,
          conditionalPanel(cond("logistic"), fluidRow(
            column(6, sl("L50", "L50: 50% retained (cm)", lm)),
            column(6, sl("L95", "L95: 95% retained (cm)", 1.2 * lm)))),
          conditionalPanel(cond("normal"), fluidRow(
            column(6, sl("Lpeak", "Best-caught length (cm)", lm)),
            column(6, sliderInput(id("sd"), "Spread (cm)", min = 0.5, max = Lmax / 2, value = v("sd", lm / 3), step = 0.5)))),
          conditionalPanel(cond("dnormal"), fluidRow(
            column(4, sliderInput(id("Lpeak_d"), "Best-caught length (cm)", min = 0, max = Lmax,
                                  value = v("Lpeak", lm), step = 0.5)),
            column(4, sliderInput(id("sdl"), "Spread, small side (cm)", min = 0.5, max = Lmax / 2, value = v("sd_left", lm / 3), step = 0.5)),
            column(4, sliderInput(id("sdr"), "Spread, large side (cm)", min = 0.5, max = Lmax, value = v("sd_right", lm / 2), step = 0.5)))))))
    })
  })
  outputOptions(output, "sel_rows", suspendWhenHidden = FALSE)

  sel_now <- reactive({
    x <- collected()
    validate(need(!inherits(x, "error"), "Check the values."))
    adapt_species(x, r$fishery(), r$one_gear())
  })

  output$sp_sel <- renderPlot({
    s <- sel_now()
    L <- seq(0, s$lh$linf * 1.2, length.out = 300)
    d <- do.call(rbind, lapply(names(s$sel), function(f) data.frame(L = L, sel = sel_len(s$sel[[f]], L), fleet = f)))
    m <- data.frame(L = L, mat = 1 / (1 + exp(-log(19) * (L - s$lh$lm) / (0.15 * s$lh$lm))))
    ggplot(d, aes(L, sel, colour = fleet)) + geom_line(linewidth = 1.1) +
      scale_colour_manual(values = app_cols(names(s$sel))) +
      geom_line(data = m, aes(L, mat), inherit.aes = FALSE, linetype = 2, colour = "grey30") +
      geom_vline(xintercept = s$lh$linf, linetype = 3, colour = "grey50") +
      labs(x = "Length (cm)", y = "Selectivity / proportion mature", colour = NULL,
           caption = "Dotted: Linf") + theme_bw(base_size = 14) + theme(legend.position = "bottom")
  })

  eq_plot <- eventReactive(input$sp_eq_go, {
    s <- sel_now()
    withProgress(message = "Calculating equilibrium curves.", value = 0.5,
                 plot_eq_selectivity(s, s$sel, fmax = 4))
  })
  output$sp_eq <- renderPlot({
    validate(need(input$sp_eq_go > 0, "Press Draw (takes a few seconds)."))
    eq_plot() + theme_bw(base_size = 14)
  })

  # ---- exploitation ----
  output$rq_rows <- renderUI({
    rv$applied
    s <- cur_adapted()
    fl <- r$fleets()
    if (r$one_gear()) return(NULL)
    lapply(seq_along(fl), function(i) fluidRow(
      column(2, div(class = "fleet-name", fl[i])),
      column(5, sliderInput(paste0("rq_", fid(fl[i])), lab("Relative catchability", "relq"), min = 0, max = 2,
                            value = round(s$relq[[fl[i]]], 2), step = 0.05))))
  })
  outputOptions(output, "rq_rows", suspendWhenHidden = FALSE)
  output$rq_onegear <- renderUI({
    if (r$one_gear()) div(class = "infobox", "One-gear fishery: catchability is set so that effort equals F/Fmsy. Set F/Fmsy with the effort trend in Fishery.")
  })
}
