# Step 6: Export (Datasets, Settings). Writes the chosen species and
# iterations with export_data() and offers them as one zip file.

export_ui <- function() {
  tagList(
    uiOutput("ex_empty"),
    tabsetPanel(
      id = "ex_tabs",
      tabPanel("Datasets",
        box(width = 12, title = "Choose what to export",
          fluidRow(
            column(6, uiOutput("ex_sp_ui")),
            column(6, uiOutput("ex_it_ui"))),
          fluidRow(
            column(4, radioButtons("ex_give", "Give trainees",
              choices = c("Data and life history" = "data_lh", "Data only" = "data"))),
            column(4, checkboxGroupInput("ex_format", "File format",
              choices = c("csv (one folder per dataset)" = "csv", "Excel (one file per dataset)" = "xlsx"),
              selected = "csv")),
            column(4, uiOutput("ex_key_ui"))),
          downloadButton("ex_download", "Download datasets (.zip)", class = "btn-primary"),
          p(class = "note", "One folder or file per dataset (iteration) with catch, effort, cpue, survey, lfd (len_mm, n) and bio (len_mm, age, wt_g, mature), plus README.txt and settings.yaml, which re-creates the data.")),
        box(width = 12, title = "What is in the export", uiOutput("ex_preview"))),
      tabPanel("Settings",
        box(width = 12, title = "Save the settings",
          p("The fishery, species and sampling settings, with the seed, as one YAML file. With the same FLOSS version, they re-create the same operating model and data."),
          downloadButton("ex_settings", "Download settings (.yaml)"),
          tags$hr(),
          p(tags$strong("In R:")),
          tags$pre("library(FLOSS)\nst  <- yaml::read_yaml(\"floss_settings.yaml\")\n# or rebuild step by step with fishery(), species(), build_oms(), oem_settings(), sample_data()")))
    )
  )
}

export_server <- function(input, output, session, rv, r) {

  output$ex_empty <- renderUI({
    if (is.null(rv$data)) empty_msg("Sample data first (Sampling step and the Sample data button).")
    else if (!r$data_current()) div(class = "warnbox", "Settings changed since the data were sampled. The export uses the last sampled data.")
  })

  output$ex_sp_ui <- renderUI({
    req(rv$data)
    sp <- unique(rv$data$catch$species)
    checkboxGroupInput("ex_sp", "Species", choices = sp, selected = sp)
  })
  output$ex_it_ui <- renderUI({
    req(rv$data)
    its <- rv$data$settings$iters
    selectizeInput("ex_it", "Datasets (iterations)", choices = its, selected = its, multiple = TRUE,
                   options = list(plugins = list("remove_button")))
  })
  output$ex_key_ui <- renderUI({
    if (r$show_truth()) checkboxInput("ex_key", "Include the answer key (answer_key/ folder)", TRUE)
    else div(class = "note", icon("eye-slash"), " The answer key can be added when Hide truth is unticked.")
  })

  sub <- reactive({
    req(rv$data, length(input$ex_sp) > 0, length(input$ex_it) > 0)
    subset_data(rv$data, input$ex_sp, as.integer(input$ex_it))
  })

  output$ex_preview <- renderUI({
    d <- sub()
    tagList(
      tags$ul(lapply(unique(d$catch$species), function(s) tags$li(tags$strong(s), ": ", data_counts(d, s)))),
      p(class = "note", paste0(length(d$settings$iters), " dataset(s): iterations ", paste(d$settings$iters, collapse = ", "),
                               ". Answer key: ", if (r$show_truth() && isTRUE(input$ex_key)) "included" else "not included", ".")))
  })

  output$ex_download <- downloadHandler(
    filename = function() paste0("floss_data_", format(Sys.Date(), "%Y%m%d"), ".zip"),
    content = function(file) {
      tmp <- file.path(tempfile("floss"), "floss_data")
      fmt <- input$ex_format %||% "csv"
      export_data(sub(), tmp, give = input$ex_give %||% "data_lh", format = fmt,
                  answer_key = r$show_truth() && isTRUE(input$ex_key), zip = TRUE)
      file.copy(paste0(tmp, ".zip"), file)
    })

  output$ex_settings <- downloadHandler(
    filename = function() "floss_settings.yaml",
    content = function(file) {
      yaml::write_yaml(list(
        floss_version = as.character(utils::packageVersion("FLOSS")),
        seed = input$fish_seed, nsim = input$nsim, datasets = input$samp_n,
        fishery = unclass_list(r$fishery()),
        species = lapply(r$species_used(), unclass_list),
        oem = unclass_list(r$oem())), file)
    })
}

# drop classes recursively (for YAML)
unclass_list <- function(x) if (is.list(x)) lapply(unclass(x), unclass_list) else x

# Keep some species and iterations of a floss_data
subset_data <- function(d, species, iters) {
  for (s in c("catch", "cpue", "survey", "lfd", "bio", "truth", "refpts", "lhpars")) {
    x <- d[[s]]
    if (is.null(x)) next
    if (!is.null(x$species)) x <- x[x$species %in% species, , drop = FALSE]
    if (!is.null(x$iter)) x <- x[x$iter %in% iters, , drop = FALSE]
    d[[s]] <- x
  }
  if (!is.null(d$effort)) d$effort <- d$effort[d$effort$iter %in% iters, , drop = FALSE]
  d$settings$iters <- iters
  keep <- vapply(d$settings$species, function(s) s$species %in% species, TRUE)
  d$settings$species <- d$settings$species[keep]
  d
}
