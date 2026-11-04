# Step 7: True status (Status, Trajectories, Kobe plot, Downloads). In the menu
# only when Hide truth is unticked: the trainer's view of the answer key.

truth_ui <- function() {
  tagList(
    uiOutput("tr_hidden"),
    tabsetPanel(
      id = "truth_tabs",
      tabPanel("Status",
        box(width = 12, title = "True status in the last year",
          uiOutput("ts_boxes"),
          tableOutput("ts_summary"),
          p(class = "note", "Median over all iterations of the operating model. Healthy: spawning biomass at or above Bmsy and F at or below Fmsy.")),
        box(width = 12, title = "Each sampled dataset",
          tableOutput("ts_table"),
          p(class = "note", "The status that an assessment of each exported dataset should find."))),
      tabPanel("Trajectories",
        box(width = 12, title = "True B/Bmsy, F/Fmsy and catch",
          plotOutput("ts_plot", height = 480),
          p(class = "note", "One line per iteration chosen at the top."))),
      tabPanel("Kobe plot",
        box(width = 12, title = "True status: Kobe plot",
          radioButtons("ts_kobe_sp", NULL, inline = TRUE,
                       choices = c("All species" = "all", "Only the species shown at the top" = "shown")),
          plotOutput("ts_kobe", height = 520),
          p(class = "note", "The iterations chosen at the top, coloured by species. One iteration: the whole history as a path from the first year (open circle) to the last (filled). Use the arrows at the top to step through iterations, 1-5 or All for many."))),
      tabPanel("Downloads",
        box(width = 12, title = "True time series and reference points (for trainers)",
          downloadButton("ts_download", "Download truth (.zip)"),
          p(class = "note", "truth_timeseries.csv (SSB, recruitment, Fbar, catch, B/Bmsy, F/Fmsy by year and iteration), refpts.csv and status.csv, for all iterations of the operating model.")))
    )
  )
}

truth_server <- function(input, output, session, rv, r) {

  output$tr_hidden <- renderUI({
    if (!r$show_truth()) return(hidden_msg())
    if (is.null(rv$oms)) empty_msg("Build the operating model first.")
  })

  oms <- reactive({ req(r$show_truth(), rv$oms); rv$oms })

  output$ts_boxes <- renderUI({
    st <- status_table(oms(), summary = TRUE)
    fluidRow(lapply(seq_len(nrow(st)), function(i) column(3,
      stat_box(st$species[i], HTML(paste0("B/Bmsy ", sprintf("%.2f", st$B_Bmsy[i]),
                                          "<br>F/Fmsy ", sprintf("%.2f", st$F_Fmsy[i])))))))
  })
  output$ts_summary <- renderTable(int_year(status_table(oms(), summary = TRUE)), digits = 2)

  output$ts_table <- renderTable({
    req(rv$data)
    st <- status_table(rv$data)
    int_year(st[st$species == r$species_used()[[r$shown()]]$species, ])
  }, digits = 2)

  output$ts_plot <- renderPlot(plot_truth(oms(), species = r$species_used()[[r$shown()]]$species,
                                          iters = r$iters_om()) + theme_bw(base_size = 14))
  output$ts_kobe <- renderPlot(plot_kobe(oms(), iters = r$iters_om(),
                                         species = if (identical(input$ts_kobe_sp, "shown"))
                                           r$species_used()[[r$shown()]]$species) + theme_bw(base_size = 14) +
                                 theme(legend.position = "bottom", legend.text = element_text(face = "italic")))

  output$ts_download <- downloadHandler(
    filename = function() paste0("floss_truth_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      d <- file.path(tempfile("truth"), "floss_truth")
      export_truth(oms(), d)
      old <- setwd(dirname(d)); on.exit(setwd(old))
      utils::zip(file, basename(d), flags = "-rq")
    })
}

# whole numbers for year and iteration columns in tables
int_year <- function(d) {
  for (n in intersect(c("year", "iter"), names(d))) d[[n]] <- as.integer(d[[n]])
  if ("p_healthy" %in% names(d)) names(d)[names(d) == "p_healthy"] <- "share healthy"
  d
}
