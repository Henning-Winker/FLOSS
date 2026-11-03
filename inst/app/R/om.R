# Step 3: Operating model (Settings, Stock trends, Fleets, Equilibrium,
# Spawners vs vulnerable, Stock and recruitment)

om_ui <- function() {
  tagList(
    page_species_ui("om_sp_ui"),
    uiOutput("om_empty"),
    tabsetPanel(
      id = "om_tabs",
      tabPanel("Settings",
        box(width = 12, title = "Operating model settings",
          fluidRow(
            column(4, sliderInput("nsim", lab("Iterations", "nsim"), min = 5, max = 100, value = 20, step = 5)),
            column(4, br(), checkboxInput("ar1", lab("Autocorrelated recruitment", "ar1"), FALSE))),
          p(class = "note", "Each iteration is one possible history of the stock with the same fishery and biology but different recruitment and effort noise. Each can become one dataset in Sampling."),
          div(class = "infobox", tags$strong("How the stock is built. "),
            "Life history gives growth, maturity and natural mortality at age. Each species is scaled so that its equilibrium maximum sustainable yield equals the target MSY, under the combined selectivity of all fleets. Fishing mortality by fleet = relative catchability x Fmsy x efficiency x effort. The stock is then projected year by year with random recruitment.")),
        box(width = 12, title = "Species in the operating model",
          tableOutput("om_table"),
          p(class = "note", "Reference points of the realised fleet mix. Shown only when Hide truth is unticked."))),
      tabPanel("Stock trends",
        box(width = 12, title = "True stock status: all species",
          uiOutput("om_hide1"),
          plotOutput("om_all", height = 520),
          p(class = "note", "Spawning biomass relative to Bmsy, fishing mortality relative to Fmsy, catch relative to MSY, and recruitment. Median and 90% interval over iterations."))),
      tabPanel("Fleets",
        box(width = 12, title = "Fishing mortality and catch by fleet",
          radioButtons("om_fleet_what", NULL, inline = TRUE,
                       choices = c("Fishing mortality (relative to Fmsy)" = "fapic", "Catch (t)" = "catch")),
          uiOutput("om_hide2"),
          plotOutput("om_fleets", height = 380),
          p(class = "note", "For the species shown at the top. Fleets with dome-shaped selectivity can take much catch at moderate F."))),
      tabPanel("Equilibrium",
        box(width = 12, title = "Production curve and the observed path",
          fluidRow(column(4, checkboxInput("om_eq_rel", "Relative to MSY reference points", FALSE)),
                   column(4, checkboxInput("om_eq_obs", "Show the simulated history (when Hide truth is unticked)", TRUE))),
          plotOutput("om_eq", height = 460),
          p(class = "note", "Equilibrium yield, spawning biomass and recruitment against fishing mortality and spawning biomass: what the stock would settle at if fished at each level for a long time. Points: the simulated years (first iteration shown at the top). Off-curve points mean the stock is not at equilibrium."))),
      tabPanel("Spawners vs vulnerable",
        box(width = 12, title = "Who spawns and who is caught",
          plotOutput("om_ssbvb", height = 380),
          div(class = "predictbox", tags$strong("Predict. "),
            "If a fleet's selectivity is dome-shaped, will its CPUE fall faster or slower than spawning biomass when fishing increases? Check with the right panel."),
          p(class = "note", "Left: selectivity at length by fleet compared with maturity. Right: spawning biomass (reproductive potential) and each fleet's vulnerable biomass (what its catch and CPUE track), relative to the first year. Where they part, an index from that fleet misrepresents the spawning stock."))),
      tabPanel("Stock and recruitment",
        box(width = 12, title = "Stock-recruitment relationship",
          plotOutput("om_srr", height = 380),
          p(class = "note", "Beverton-Holt curve with steepness from the species, and the simulated recruitment of three iterations. Scatter around the curve is recruitment variability (and autocorrelation if ticked).")))
    )
  )
}

om_server <- function(input, output, session, rv, r) {

  output$om_empty <- renderUI({
    if (is.null(rv$oms)) empty_msg("Build the operating model first, using the Build operating model button in the left panel.")
    else if (!r$om_current()) div(class = "warnbox", "The fishery or species changed since the operating model was built. The plots show the last build. Press Build operating model to update.")
  })

  output$om_table <- renderTable({
    req(rv$oms)
    d <- do.call(rbind, lapply(rv$oms, function(o) {
      rp <- o$stock@refpts
      x <- data.frame(Species = o$species$species, Iterations = dims(o$stock)$iter,
                      Years = paste(range(dimnames(o$stock)$year), collapse = "-"))
      if (r$show_truth())
        x <- cbind(x, data.frame(Fmsy = round(c(rp["Fmsy"]), 3), `Bmsy (t)` = round(c(rp["Bmsy"])),
                                 `MSY (t)` = round(c(rp["MSY"])), `B0 (t)` = round(c(rp["B0"])),
                                 check.names = FALSE))
      x
    }))
    d
  })

  output$om_hide1 <- renderUI(if (!r$show_truth()) hidden_msg())
  output$om_hide2 <- renderUI(if (!r$show_truth()) hidden_msg())

  output$om_all <- renderPlot({
    req(rv$oms, r$show_truth())
    plot_om(rv$oms) + theme_bw(base_size = 14)
  })

  output$om_fleets <- renderPlot({
    req(r$show_truth())
    plot_fleets(r$om(), what = input$om_fleet_what %||% "fapic") + theme_bw(base_size = 14)
  })

  output$om_eq <- renderPlot({
    plot_eq(r$om(), iter = r$iters()[1], obs = isTRUE(input$om_eq_obs) && r$show_truth(), rel = isTRUE(input$om_eq_rel))
  })

  output$om_ssbvb <- renderPlot(plot_ssb_vb(r$om()) )

  output$om_srr <- renderPlot(plot_srr(r$om(), iters = 1:min(3, dims(r$om()$stock)$iter)) + theme_bw(base_size = 14))
}
