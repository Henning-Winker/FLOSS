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
          p(class = "note", "The iterations chosen at the top are drawn as lines over the median and 90% interval of all iterations."))),
      tabPanel("Kobe plot",
        box(width = 12, title = "Status of every iteration",
          plotOutput("ts_kobe", height = 460))),
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
                                          iters = r$iters()) + theme_bw(base_size = 14))
  output$ts_kobe <- renderPlot(plot_kobe(oms()) + theme_bw(base_size = 14))

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
