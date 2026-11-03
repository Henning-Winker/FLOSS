# The layout: left panel (steps, Hide truth, Build/Sample block) and the page
# for each step. Same structure and style as the FLicc app.

app_sidebar <- function() {
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      id = "step",
      menuItem("Guide", tabName = "guide", icon = icon("book-open")),
      menuItem("1  Fishery *", tabName = "fishery", selected = TRUE),
      menuItem("2  Species *", tabName = "species"),
      menuItem("3  Operating model *", tabName = "om"),
      menuItem("4  Treatments", tabName = "treat"),
      menuItem("5  Sampling *", tabName = "sampling"),
      menuItem("6  Export *", tabName = "export"),
      menuItemOutput("truth_item")
    ),
    div(style = "padding: 4px 16px;",
      tags$p(style = "font-size:12px; color:#8aa4af;", "* required. All other steps are optional."),
      tags$hr(style = "border-color:#4b646f; margin: 8px 0;"),
      div(class = "truth-switch",
          checkboxInput("hide_truth", "Hide truth", value = FALSE),
          tags$p(style = "font-size:12px; color:#8aa4af; margin-top:-8px;",
                 "When ticked, true stock status and the answer key are hidden (for trainees). Untick to show the True status step.")),
      tags$hr(style = "border-color:#4b646f; margin: 8px 0;"),
      tags$div(style = "font-size:12px; font-weight:600; color:#8aa4af; text-transform:uppercase;", "Before you run"),
      uiOutput("checklist"),
      actionButton("build_om", "Build operating model", class = "btn-primary", width = "100%",
                   style = "margin-top:10px; height:44px; font-weight:600;"),
      actionButton("sample", "Sample data", class = "btn-primary", width = "100%",
                   style = "margin-top:8px; height:44px; font-weight:600;"),
      tags$p(style = "margin-top:8px; font-size:13px;", textOutput("run_help", inline = TRUE)),
      actionButton("start_over", "Start over", width = "100%",
                   style = "margin-top:10px; height:40px; background:transparent; color:#7fc0e8; border:1px solid #5d7480;")
    )
  )
}

app_css <- "
  .box { border-top-color: #3c8dbc; }
  .content-wrapper { background: #ecf0f5; }
  .skin-blue .main-header .logo { font-weight: 700; }
  .note { font-size: 13px; color: #555; margin-top: 6px; }
  .warnbox { border: 1px solid #f5a623; background: #fff4e0; color: #5c3b00; padding: 10px 12px; border-radius: 3px; margin: 8px 0; font-size: 14px; }
  .infobox { border: 1px solid #9cc3dc; background: #eef6fb; color: #22404f; padding: 10px 12px; border-radius: 3px; margin: 8px 0; font-size: 14px; }
  .predictbox { border: 1px solid #b7a6dd; background: #f4f0fb; color: #2f2350; padding: 10px 12px; border-radius: 3px; margin: 8px 0; font-size: 14px; }
  .flagbox { border: 1px solid #e3a0a0; background: #fcf0f0; color: #5a1f1f; padding: 10px 12px; border-radius: 3px; margin: 8px 0; font-size: 14px; }
  .statbox { background:#fff; border-radius:3px; box-shadow:0 1px 1px rgba(0,0,0,0.1); padding:10px 14px; margin-bottom:10px; }
  .statbox .lab { font-size:13px; color:#555; }
  .statbox .val { font-size:26px; font-weight:700; color:#333; }
  .formatbox { border: 1px dashed #b9c3c9; background: #fafbfc; padding: 10px 12px; font-size: 13px; margin-top: 10px; }
  .tip { color: #3c8dbc; cursor: help; margin-left: 4px; font-size: 12px; }
  .fleet-row { border-bottom: 1px solid #eef1f3; padding-top: 6px; margin-bottom: 4px; }
  .fleet-name { font-weight: 700; padding-top: 28px; }
  .truth-switch .checkbox label { font-weight: 600; color: #fff; }
  .irs--shiny .irs-bar { background: #3c8dbc; border-color: #3c8dbc; }
  .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background: #3c8dbc; }
  .guide h2 { font-size: 20px; margin-top: 22px; }
  .guide table { border-collapse: collapse; margin: 8px 0; }
  .guide th, .guide td { border: 1px solid #d8dee2; padding: 4px 8px; font-size: 13px; }
  .guide pre { background: #f6f8fa; border: none; }
"

app_body <- function() {
  dashboardBody(
    tags$head(tags$style(HTML(app_css))),
    div(style = "display:flex; flex-wrap:wrap; align-items:center; justify-content:flex-end; gap:14px; margin-bottom:6px;",
        uiOutput("species_selector"), uiOutput("iter_selector"), uiOutput("status_pill")),
    tabItems(
      tabItem("guide", guide_ui()),
      tabItem("fishery", fishery_ui()),
      tabItem("species", species_ui()),
      tabItem("om", om_ui()),
      tabItem("treat", treat_ui()),
      tabItem("sampling", sampling_ui()),
      tabItem("export", export_ui()),
      tabItem("truth", truth_ui())
    )
  )
}
