# Guide: the user guide (inst/guide/how-it-works.md, the same text as the
# vignette "How FLOSS works"), a teaching flow and the settings explained.

guide_sections <- function() {
  f <- system.file("guide", "how-it-works.md", package = "FLOSS")
  if (!nzchar(f)) return(list())
  x <- readLines(f, warn = FALSE)
  starts <- grep("^## ", x)
  ends <- c(starts[-1] - 1, length(x))
  secs <- Map(function(a, b) paste(x[a:b], collapse = "\n"), starts, ends)
  names(secs) <- sub("^## [0-9]+\\. *", "", x[starts])
  secs
}

md_html <- function(txt) HTML(commonmark::markdown_html(txt, extensions = TRUE))

guide_ui <- function() {
  s <- guide_sections()
  grp <- function(nm) div(class = "guide", md_html(paste(unlist(s[nm]), collapse = "\n\n")))
  tabsetPanel(
    id = "guide_tabs",
    tabPanel("Start here",
      box(width = 12, title = "What FLOSS is for",
        p("FLOSS simulates fisheries whose truth is known, and the data an observer would collect from them.",
          "It is built for teaching population dynamics: predict how a stock responds to a change, see the response in the true stock and in the data, and learn which patterns in data are warning flags rather than real dynamics."),
        div(class = "infobox", tags$strong("The steps. "),
          "Work down the left panel: ", tags$strong("1 Fishery"), " (fleets and effort), ",
          tags$strong("2 Species"), " (life history and selectivity), ",
          tags$strong("3 Operating model"), " (press Build operating model), ",
          tags$strong("4 Treatments"), " (optional: predict, then reveal), ",
          tags$strong("5 Sampling"), " (press Sample data) and ", tags$strong("6 Export"), "."),
        div(class = "predictbox", tags$strong("Predict first. "),
          "Before you look at a result, write down what you expect: up, down or no change, and when. The Treatments step asks for your prediction before it shows the answer."),
        div(class = "flagbox", tags$strong("Warning flags. "),
          "Some patterns in data cannot come from real stock dynamics: an index that keeps rising while catches fall with steady effort, mean length jumping in one year, fish older than the plus group. FLOSS can plant such problems on purpose (Sampling step), and the answer key lists them."),
        div(class = "formatbox", tags$strong("Hide truth. "),
          "Unticked by default, so the truth is shown: the True status step, the true lines on the data plots and the answer key in the export. Tick it before showing the app to trainees, so they see the data only.")),
      box(width = 12, title = "Scope",
        p("FLOSS covers dynamics and data. Fitting assessment models to the data follows afterwards, e.g. length-based methods in FLicc and production models in JABBA, using the exported datasets and the answer key to judge the results."),
        p(class = "note", "Structure inspired by Haddon, M. (2021) Using R for Modelling and Quantitative Methods in Fisheries. CRC Press."))),
    tabPanel("How it works", box(width = 12, grp(names(s)[1:5]))),
    tabPanel("Sampling and size data", box(width = 12, grp(names(s)[6:7]))),
    tabPanel("Warning flags", box(width = 12, grp(names(s)[8]))),
    tabPanel("Truth and export", box(width = 12, grp(names(s)[9]))),
    tabPanel("Settings explained",
      box(width = 12, title = "What each setting means",
        tableOutput("guide_tips"),
        p(class = "note", "The same texts appear when you hover over the (i) symbols.")))
  )
}

guide_tips_table <- function() {
  t <- floss_tip()
  data.frame(Setting = names(t), Meaning = unlist(t), row.names = NULL)
}
