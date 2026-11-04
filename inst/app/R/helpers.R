# Small helper functions and the shared reactive objects used across the app.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Use FLOSS's own functions even when another package attached later has a
# function of the same name (e.g. JABBA3::status_table(), FLicc::plot_lfd()):
# copy the FLOSS exports into the app's environment.
for (.f in getNamespaceExports("FLOSS")) assign(.f, getExportedValue("FLOSS", .f))
rm(.f)

# ---- small UI pieces -----------------------------------------------------------

# An (i) symbol with the explanation from floss_tip() shown on hover
tip <- function(id) tags$span(class = "tip", title = floss_tip(id), icon("circle-info"))

# A label with a tip
lab <- function(text, id = NULL) if (is.null(id)) text else tags$span(text, tip(id))

# A simple number box
stat_box <- function(label, value) {
  div(class = "statbox", div(class = "lab", label), div(class = "val", value))
}

# Message shown on a page that needs something done first
empty_msg <- function(text) div(class = "infobox", text)

hidden_msg <- function() {
  div(class = "infobox", icon("eye-slash"), " Hidden while ", tags$strong("Hide truth"),
      " is ticked (left panel). Untick it to show the true stock.")
}

# ---- fishery presets -------------------------------------------------------------

fishery_presets <- list(
  reef = list(
    label = "Reef: traps, gillnet, linefish",
    fleets = c("traps", "gillnet", "linefish"),
    trends = list(
      traps    = list(e0 = 0.05, fhi = 1.0, flo = 0.8, breaks = c(0.4, 0.6)),
      gillnet  = list(e0 = 0.05, fhi = 0.5, flo = 1.5, breaks = c(0.5, 0.8)),
      linefish = list(e0 = 0.05, fhi = 2.5, flo = 0.9, breaks = c(0.6, 0.8))),
    tail = c(traps = 5, gillnet = 10, linefish = 15),
    days = c(traps = 4000, gillnet = 1500, linefish = 3000),
    rho = c(traps = 0.3, gillnet = 0.3, linefish = 0.6),
    ess_age = c(traps = 30, gillnet = 30, linefish = 50),
    ess_len = c(traps = 300, gillnet = 300, linefish = 500),
    bio_n = c(traps = 100, gillnet = 100, linefish = 200)),
  one = list(
    label = "One gear: handline (effort = F/Fmsy)",
    fleets = "handline",
    trends = list(handline = list(e0 = 0.05, fhi = 1.5, flo = 1.0, breaks = c(0.4, 0.6))),
    tail = c(handline = 15), rho = c(handline = 0.3), days = c(handline = 2000),
    ess_age = c(handline = 50), ess_len = c(handline = 500),
    bio_n = c(handline = 200))
)

# A new fishery: the gears the user added, with default settings for each
custom_preset <- function(gears) {
  n <- length(gears)
  rep_named <- function(v) stats::setNames(rep(v, n), gears)
  list(label = "New fishery: add your own gears", fleets = gears,
       trends = stats::setNames(rep(list(list(e0 = 0.05, fhi = 1, flo = 0.8, breaks = c(0.4, 0.6))), n), gears),
       tail = rep_named(10), rho = rep_named(0.3), ess_age = rep_named(30), ess_len = rep_named(300),
       bio_n = rep_named(100), days = rep_named(1000))
}

preset <- function(input, rv = NULL) {
  if (identical(input$fish_preset, "custom")) return(custom_preset(rv$gears %||% character(0)))
  fishery_presets[[input$fish_preset %||% "reef"]]
}

# Input id suffix for a gear: inputs are keyed by gear name, not position, so
# adding or removing a gear keeps the settings of the others
fid <- function(x) gsub("[^A-Za-z0-9]", "_", x)

# Current value of an input if it exists, else the default (keeps settings
# when rows are drawn again)
keepv <- function(input, id, default) {
  v <- isolate(input[[id]])
  if (is.null(v) || (length(default) == 2 && length(v) != 2)) default else v
}

# The species shown first: the two-spot red snapper when it is in the model
default_species <- function(ids) if ("lutjanus_bohar" %in% ids) "lutjanus_bohar" else ids[1]

# Species dropdown for a page, kept in step with the "Showing" selector
page_species_ui <- function(id) {
  div(style = "max-width:360px; margin-bottom:4px;", uiOutput(id))
}

# ---- species -----------------------------------------------------------------------

# The species library (four examples plus the user's saved species)
initial_species <- function() {
  lib <- species_library()
  sp <- lapply(lib$id, load_species)
  names(sp) <- lib$id
  sp
}

species_label <- function(sp) {
  paste0(sp$species, if (nzchar(sp$common_name %||% "")) paste0(" (", sp$common_name, ")"))
}

# Fit a species to the fleets of the current fishery: missing fleets get a
# logistic selectivity at the length at maturity; with the one-gear preset,
# catchability is scaled so that effort equals F/Fmsy.
adapt_species <- function(sp, fishery, one_gear = FALSE) {
  fl <- fishery$fleets
  if (one_gear) return(one_gear_species(sp, fl[1], sel = sp$sel[[fl[1]]]))
  sel <- lapply(fl, function(f) sp$sel[[f]] %||%
    selectivity("logistic", L50 = sp$lh$lm, L95 = 1.2 * sp$lh$lm))
  relq <- vapply(fl, function(f) sp$relq[[f]] %||% 0.5, 0)
  names(sel) <- names(relq) <- fl
  species(sp$species, sp$common_name, lh = sp$lh, msy = sp$msy, sel = sel,
          relq = relq, productivity = sp$productivity, id = sp$id)
}

# Selectivity at length (cm), for the previews. Shapes as in FLRef::sel_la()
sel_len <- function(s, L) {
  switch(s$type,
    logistic = 1 / (1 + exp(-log(19) * (L - s$L50) / (s$L95 - s$L50))),
    normal   = exp(-(L - s$Lpeak)^2 / (2 * s$sd^2)),
    dnormal  = ifelse(L < s$Lpeak, exp(-(L - s$Lpeak)^2 / (2 * s$sd_left^2)),
                                   exp(-(L - s$Lpeak)^2 / (2 * s$sd_right^2))))
}

# Fleet colours (the FLOSS palette): fleets in fishery order, then the survey
app_cols <- function(fleets) FLOSS:::fleet_colours(fleets)

# Iteration colours (used where gears are in separate panels)
app_iter_cols <- function(iters) FLOSS:::iter_colours(iters)

sel_types <- c("S-shaped (logistic)" = "logistic", "Symmetric dome (normal)" = "normal",
               "Dome, two spreads (double normal)" = "dnormal")

# ---- shared reactive objects ------------------------------------------------------

make_reactives <- function(input, rv) {
  r <- list()

  r$one_gear <- reactive(identical(input$fish_preset, "one"))

  r$fleets <- reactive(preset(input, rv)$fleets)

  # The fishery, built from the Fishery inputs (preset values until edited)
  r$fishery <- reactive({
    p <- preset(input, rv)
    fl <- p$fleets
    if (!length(fl)) return(structure(list(msg = "Add at least one gear (New fishery)."), class = "floss_error"))
    yrs <- input$fish_years %||% c(1969, 2025)
    y0 <- yrs[1]; y1 <- yrs[2]
    trends <- lapply(seq_along(fl), function(i) {
      d <- p$trends[[fl[i]]]
      brk <- input[[paste0("fl_brk_", fid(fl[i]))]]
      brk <- if (length(brk) == 2) (brk - y0) / (y1 - y0) else d$breaks
      brk <- pmin(pmax(brk, 0.01), 0.99)
      if (brk[2] <= brk[1]) brk[2] <- min(0.99, brk[1] + 0.01)
      effort_trend(e0 = input[[paste0("fl_e0_", fid(fl[i]))]] %||% d$e0,
                   fhi = input[[paste0("fl_fhi_", fid(fl[i]))]] %||% d$fhi,
                   flo = input[[paste0("fl_flo_", fid(fl[i]))]] %||% d$flo, breaks = brk)
    })
    names(trends) <- fl
    creep <- vapply(seq_along(fl), function(i) (input[[paste0("fl_creep_", fid(fl[i]))]] %||% 0) / 100, 0)
    cstart <- vapply(seq_along(fl), function(i) as.integer(input[[paste0("fl_cstart_", fid(fl[i]))]] %||% y0), 0L)
    tryCatch(
      fishery(fleets = fl, years = seq(y0, y1), trends = trends,
              sigmaE = input$fish_sigmaE %||% 0.15, creep = creep,
              creep_start = cstart, name = input$fish_preset %||% "reef",
              seed = input$fish_seed %||% 123),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "floss_error"))
  })

  r$fishery_ok <- reactive(inherits(r$fishery(), "floss_fishery"))

  # The species in the operating model, adapted to the fleets
  r$species_used <- reactive({
    req(r$fishery_ok())
    ids <- intersect(input$sp_use %||% names(rv$sp), names(rv$sp))
    out <- lapply(ids, function(id) tryCatch(
      adapt_species(rv$sp[[id]], r$fishery(), r$one_gear()),
      error = function(e) NULL))
    names(out) <- ids
    out[!vapply(out, is.null, TRUE)]
  })

  # The species shown on the pages (top selector)
  r$shown <- reactive({
    ids <- names(r$species_used())
    req(length(ids))
    if (!is.null(input$shown) && input$shown %in% ids) input$shown else default_species(ids)
  })

  r$om <- reactive({
    req(rv$oms)
    validate(need(r$shown() %in% names(rv$oms),
      "This species is not in the operating model yet. Press Build operating model."))
    rv$oms[[r$shown()]]
  })

  r$om_snap <- reactive({
    list(fishery = r$fishery(), species = r$species_used(),
         nsim = input$nsim %||% 20, ar1 = isTRUE(input$ar1))
  })

  r$om_current <- reactive(!is.null(rv$oms) && identical(r$om_snap(), rv$om_snap))

  # Sampling settings, built from the Sampling inputs
  r$oem <- reactive({
    req(r$fishery_ok())
    p <- preset(input, rv)
    fl <- p$fleets
    f <- r$fishery()
    pf <- function(id, def, scale = 1) {
      v <- vapply(seq_along(fl), function(i)
        as.numeric(input[[paste0(id, "_", fid(fl[i]))]] %||% def[[min(i, length(def))]]) * scale, 0)
      stats::setNames(v, fl)
    }
    last <- f$years[2]
    years <- lapply(seq_along(fl), function(i)
      as.integer(input[[paste0("sm_years_", fid(fl[i]))]] %||% c(last - p$tail[[i]] + 1, last)))
    names(years) <- fl
    sv <- NULL
    if (isTRUE(input$sv_on)) {
      sv <- survey_settings(
        years = as.integer(input$sv_years %||% c(last - 24, last)),
        timing = input$sv_timing %||% 0, cv = input$sv_cv %||% 0.25,
        rho = input$sv_rho %||% 0,
        sel = survey_sel_input(input, rv),
        lfd = TRUE, ess_age = input$sv_ess_age %||% 50, ess_len = input$sv_ess_len %||% 500,
        q_step = if (isTRUE(input$sv_qstep_on)) c(input$sv_qstep_year, input$sv_qstep_mult))
    }
    qstep <- NULL
    if (isTRUE(input$cp_qstep_on) && !is.null(input$cp_qstep_fleet)) {
      qstep <- stats::setNames(list(c(input$cp_qstep_year, input$cp_qstep_mult)), input$cp_qstep_fleet)
    }
    seasons <- tryCatch(app_seasons(input, fl), error = function(e) e)
    if (inherits(seasons, "error"))
      return(structure(list(msg = conditionMessage(seasons)), class = "floss_error"))
    tryCatch(oem_settings(f, seasons = seasons,
      years = years,
      catch_cv = pf("ct_cv", 0.1), catch_bias = pf("ct_bias", 0, 1 / 100),
      catch_bias_from = if (isTRUE(input$ct_from_on)) input$ct_from,
      effort_cv = pf("ef_cv", 0.1), effort_days = pf("ef_days", p$days),
      cpue_cv = pf("cp_cv", 0.2), cpue_rho = pf("cp_rho", p$rho), cpue_beta = pf("cp_beta", 1),
      cpue_creep = isTRUE(input$cp_creep %||% TRUE), cpue_q_step = qstep,
      lfd_ess_age = pf("lf_ess_age", p$ess_age), lfd_ess_len = pf("lf_ess_len", p$ess_len),
      lfd_bin = if (isTRUE((input$lf_bin %||% 0) > 0)) input$lf_bin,
      bio_n = pf("bi_n", p$bio_n), bio_p_aged = input$bi_p_aged %||% 0.5,
      bio_age_scheme = input$bi_scheme %||% "random", bio_n_per_bin = input$bi_n_per_bin %||% 5,
      age_cv = input$bi_age_cv %||% 0.1, bio_age_bias = (input$bi_age_bias %||% -5) / 100,
      wt_cv = input$bi_wt_cv %||% 0.1, survey = sv),
      error = function(e) structure(list(msg = conditionMessage(e)), class = "floss_error"))
  })

  r$data_snap <- reactive(list(om = rv$om_snap, oem = r$oem(), n = input$samp_n %||% 5))

  r$data_current <- reactive(!is.null(rv$data) && r$om_current() &&
                               identical(r$data_snap(), rv$data_snap))

  # Iterations shown on the pages (top selector)
  # Iterations chosen at the top: any of the operating model's iterations
  # (truth, operating model pages); data pages use those that were sampled
  r$nsim <- reactive(if (is.null(rv$oms)) 0L else dims(rv$oms[[1]]$stock)$iter)
  r$iters_om <- reactive({
    v <- as.integer(input$iters_shown %||% 1)
    v <- v[v >= 1 & v <= max(1, r$nsim())]
    if (!length(v)) 1L else v
  })
  r$iters <- reactive({
    its <- rv$data$settings$iters %||% 1
    v <- r$iters_om()
    v <- v[v %in% its]
    if (!length(v)) its[1] else v
  })

  # The answer key can be shown?
  r$show_truth <- reactive(!isTRUE(input$hide_truth))

  # Can the operating model be built? (and what is missing)
  r$ready_om <- reactive({
    miss <- character(0)
    if (!r$fishery_ok()) miss <- c(miss, "fix the fishery settings")
    else if (!length(r$species_used())) miss <- c(miss, "choose at least one species")
    list(ok = length(miss) == 0, missing = miss)
  })

  r
}

# Seasonal sampling design from the Seasons inputs: the design for the chosen
# fleets, balanced for the others
app_seasons <- function(input, fl) {
  n <- as.integer(input$se_n %||% 4)
  type <- input$se_type %||% "balanced"
  sh <- function(prefix) {
    v <- vapply(seq_len(n), function(k) as.numeric(input[[paste0(prefix, k)]] %||% 1), 0)
    if (sum(v) <= 0) stop("season shares must not all be 0")
    v
  }
  nm <- if (isTRUE(input$se_nm_on))
    list(from = input$se_nm_year, mult = input$se_nm_mult %||% 1, over = input$se_nm_over %||% 0)
  pm <- input$se_pmiss %||% 0
  mk <- function(...) season_design(n, ..., n_mult = nm, p_missing = pm)
  des <- switch(type,
    balanced = mk(),
    fixed = mk(share = sh("se_share_")),
    shift = mk(share = sh("se_share_"),
               shift = list(from = input$se_shift_year, share = sh("se_share2_"),
                            over = input$se_shift_over %||% 0)),
    random = mk(share = sh("se_share_"), kappa = input$se_kappa %||% 5))
  use <- input$se_fleets %||% fl
  out <- lapply(stats::setNames(fl, fl), function(f) if (f %in% use) des else season_design(n))
  out
}

# Survey selectivity (lengths as shares of Linf, as survey_settings() takes
# them); the sliders show them in cm for the species shown
survey_sel_input <- function(input, rv) {
  v <- rv$sv
  switch(input$sv_type %||% "logistic",
    logistic = list(type = "logistic", L50 = v$L50, L95 = max(v$L95, v$L50 + 0.01)),
    normal = list(type = "normal", Lpeak = v$Lpeak, sd = v$sd),
    dnormal = list(type = "dnormal", Lpeak = v$Lpeak, sd_left = v$sd_left, sd_right = v$sd_right))
}

# Short text for one species' data, e.g. in previews
data_counts <- function(d, sp) {
  n <- function(x) if (is.null(x)) 0 else sum(x$species == sp)
  paste0(n(d$catch), " catch, ", n(d$cpue), " CPUE, ", n(d$survey), " survey, ",
         n(d$lfd), " length-bin and ", n(d$bio), " fish records")
}
