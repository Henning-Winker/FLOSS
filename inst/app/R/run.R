# Left panel (checklist, Build operating model, Sample data), the selectors and
# status at the top of the page, and the True status menu item.

run_server <- function(input, output, session, rv, r) {

  # ---- left panel ----
  output$checklist <- renderUI({
    item <- function(state, text) {
      mark <- switch(state, ok = "✓", todo = "✗", "○")
      col <- switch(state, ok = "#6fdca8", todo = "#f5a623", "#8aa4af")
      tags$div(style = "font-size:13px; margin:3px 0;",
               tags$span(style = paste0("color:", col, "; font-weight:700; margin-right:8px;"), mark), text)
    }
    nsp <- tryCatch(length(r$species_used()), error = function(e) 0)
    tagList(
      item(if (r$fishery_ok()) "ok" else "todo", "Fishery set *"),
      item(if (nsp > 0) "ok" else "todo", paste0("Species chosen (", nsp, ") *")),
      item(if (r$om_current()) "ok" else "todo", "Operating model built *"),
      item(if (!is.null(rv$treat)) "ok" else "opt", "Treatment run (optional)"),
      item(if (r$data_current()) "ok" else "todo", "Data sampled *")
    )
  })

  output$run_help <- renderText({
    rd <- r$ready_om()
    if (!rd$ok) return(paste0("To build: ", paste(rd$missing, collapse = ", "), "."))
    if (is.null(rv$oms)) return("Operating model not built yet.")
    if (!r$om_current()) return("Out of date: fishery or species changed since the last build.")
    if (is.null(rv$data)) return("Operating model up to date. Data not sampled yet.")
    if (!r$data_current()) return("Out of date: sampling settings changed since data were sampled.")
    "Operating model and data up to date."
  })

  # ---- top of the page ----
  output$species_selector <- renderUI({
    sps <- tryCatch(r$species_used(), error = function(e) list())
    req(length(sps) > 0)
    ch <- stats::setNames(names(sps), vapply(sps, `[[`, "", "species"))
    div(style = "display:flex; align-items:center; gap:8px;",
        tags$strong("Showing:"),
        div(style = "width:230px; margin-bottom:-15px;",
            selectInput("shown", NULL, choices = ch,
                        selected = isolate(input$shown) %||% default_species(names(sps)))))
  })

  # species dropdowns on the Operating model and Sampling pages: same choice
  # as "Showing" at the top; changing one changes all
  for (pid in c("om_sp", "sm_sp")) local({
    id <- pid
    output[[paste0(id, "_ui")]] <- renderUI({
      sps <- tryCatch(r$species_used(), error = function(e) list())
      req(length(sps) > 0)
      ch <- stats::setNames(names(sps), vapply(sps, `[[`, "", "species"))
      selectInput(id, "Species", choices = ch, selected = isolate(r$shown()))
    })
    observeEvent(input[[id]], {
      if (!identical(input[[id]], input$shown)) updateSelectInput(session, "shown", selected = input[[id]])
    }, ignoreInit = TRUE)
  })
  observeEvent(input$shown, {
    for (id in c("om_sp", "sm_sp"))
      if (!identical(input[[id]], input$shown)) updateSelectInput(session, id, selected = input$shown)
  })

  # Iteration picker: all iterations of the operating model (sampled ones
  # marked with *); arrows step through single iterations
  output$iter_selector <- renderUI({
    n <- r$nsim()
    req(n > 0)
    sampled <- rv$data$settings$iters %||% integer(0)
    ch <- stats::setNames(seq_len(n), paste0(seq_len(n), ifelse(seq_len(n) %in% sampled, " *", "")))
    btn <- function(id, lab, ttl) actionButton(id, lab, title = ttl,
      style = "padding:2px 8px; margin:0 1px; font-size:13px;")
    div(style = "display:flex; align-items:center; gap:6px;",
        tags$strong("Iteration:"), tip("nsim"),
        btn("iter_prev", "\u25C0", "Previous iteration"),
        div(style = "width:170px; margin-bottom:-15px;",
            selectizeInput("iters_shown", NULL, choices = ch, multiple = TRUE,
                           selected = isolate(input$iters_shown) %||% 1,
                           options = list(placeholder = "1", plugins = list("remove_button")))),
        btn("iter_next", "\u25B6", "Next iteration"),
        btn("iter_5", "1-5", "Show iterations 1 to 5"),
        btn("iter_all", "All", "Show all iterations"))
  })
  step_iter <- function(d) {
    n <- r$nsim(); cur <- max(r$iters_om())
    updateSelectizeInput(session, "iters_shown", selected = ((cur - 1 + d) %% n) + 1)
  }
  observeEvent(input$iter_prev, step_iter(-1))
  observeEvent(input$iter_next, step_iter(1))
  observeEvent(input$iter_5, updateSelectizeInput(session, "iters_shown", selected = seq_len(min(5, r$nsim()))))
  observeEvent(input$iter_all, updateSelectizeInput(session, "iters_shown", selected = seq_len(r$nsim())))

  output$status_pill <- renderUI({
    txt <- if (!r$fishery_ok()) "Fishery settings need fixing"
      else if (is.null(rv$oms)) "Operating model not built yet"
      else if (!r$om_current()) "Operating model out of date"
      else if (is.null(rv$data)) "Data not sampled yet"
      else if (!r$data_current()) "Data out of date"
      else "Up to date"
    stale <- txt %in% c("Operating model out of date", "Data out of date", "Fishery settings need fixing")
    tags$span(style = paste0("padding:4px 12px; border-radius:12px; font-size:13px; font-weight:600; border:1px solid ",
      if (stale) "#f5a623; background:#fff4e0; color:#7a4b00;" else "#9cc3dc; background:#eef6fb; color:#1f5f87;"), txt)
  })

  # True status appears in the menu only when Hide truth is unticked
  output$truth_item <- renderMenu({
    if (isTRUE(input$hide_truth)) return(NULL)
    menuItem("7  True status", tabName = "truth", icon = icon("eye"))
  })
  observeEvent(input$hide_truth, {
    if (isTRUE(input$hide_truth) && identical(input$step, "truth"))
      updateTabItems(session, "step", "fishery")
  })

  # ---- Build operating model ----
  observeEvent(input$build_om, {
    rd <- r$ready_om()
    if (!rd$ok) {
      showNotification(paste0("To build: ", paste(rd$missing, collapse = ", "), "."), type = "warning")
      return()
    }
    sps <- r$species_used()
    f <- r$fishery()
    n <- input$nsim %||% 20
    oms <- NULL
    withProgress(message = "Building the operating model.", value = 0, {
      # species already built with the same fishery and settings are reused
      old <- if (!is.null(rv$om_snap) && identical(rv$om_snap$fishery, f) &&
                 identical(rv$om_snap$nsim, n) && identical(rv$om_snap$ar1, isTRUE(input$ar1)))
        rv$oms else NULL
      eff <- effort_sim(f, n)
      oms <- tryCatch({
        out <- list()
        for (id in names(sps)) {
          incProgress(1 / length(sps), detail = sps[[id]]$species)
          same <- !is.null(old[[id]]) && identical(old[[id]]$species, sps[[id]])
          out[[id]] <- if (same) old[[id]] else
            build_om(sps[[id]], f, nsim = n, effort = eff, ar1 = isTRUE(input$ar1))
        }
        structure(out, class = "floss_oms")
      }, error = function(e) e)
    })
    if (inherits(oms, "error")) {
      showNotification(paste("The operating model could not be built:", conditionMessage(oms)),
                       type = "error", duration = 15)
      return()
    }
    rv$oms <- oms
    rv$om_snap <- r$om_snap()
    rv$data <- NULL
    rv$treat <- NULL
  })

  # ---- Sample data ----
  observeEvent(input$sample, {
    if (!r$om_current()) {
      showNotification("Build the operating model first (or again: it is out of date).", type = "warning")
      return()
    }
    oem <- r$oem()
    if (inherits(oem, "floss_error")) {
      showNotification(paste("Check the sampling settings:", oem$msg), type = "error", duration = 12)
      return()
    }
    n <- min(input$samp_n %||% 5, input$nsim %||% 20)
    d <- NULL
    withProgress(message = "Sampling data.", value = 0.5, {
      d <- tryCatch(sample_data(rv$oms, oem, iters = seq_len(n), seed = input$fish_seed %||% 123),
                    error = function(e) e)
    })
    if (inherits(d, "error")) {
      showNotification(paste("Data could not be sampled:", conditionMessage(d)), type = "error", duration = 15)
      return()
    }
    rv$data <- d
    rv$data_snap <- r$data_snap()
    if (!identical(input$step, "sampling")) updateTabItems(session, "step", "sampling")
  })

  observeEvent(input$start_over, session$reload())
}
