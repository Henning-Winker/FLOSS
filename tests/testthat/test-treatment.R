test_that("treatments share random numbers with the baseline", {
  skip_on_cran()
  tr <- list(
    pulse = treatment("pulse", recruitment = list(type = "pulse", years = 2010, mult = 5)),
    cut   = treatment("cut", effort = list(fleet = "linefish", from = 2015, mult = 0.5)))
  rt <- run_treatment("lutjanus_bohar", reef_fishery(), tr, nsim = 3,
                      oem = reef_oem(), iters = 1)
  b <- rt$oms$baseline$lutjanus_bohar$stock
  p <- rt$oms$pulse$lutjanus_bohar$stock
  c <- rt$oms$cut$lutjanus_bohar$stock
  # identical before the treatment starts
  expect_equal(c(ssb(b)[, ac(1971:2009)]), c(ssb(p)[, ac(1971:2009)]))
  expect_equal(c(catch(b)[, ac(1971:2014)]), c(catch(c)[, ac(1971:2014)]))
  # the treatment has the intended effect
  expect_equal(c(rec(p)[, "2010"] / rec(b)[, "2010"]), rep(5, 3), tolerance = 1e-6)
  expect_true(all(c(fbar(c)[, "2020"]) < c(fbar(b)[, "2020"])))
  # sampled data identical before the change, too
  d0 <- rt$data$baseline$catch
  d1 <- rt$data$cut$catch
  k <- d0$year < 2015
  expect_equal(d0$catch_t[k], d1$catch_t[k])
})

test_that("one-gear scenario sets F/Fmsy through effort", {
  skip_on_cran()
  sp <- one_gear_species(load_species("lutjanus_bohar"))
  om <- build_om(sp, one_gear_fishery(2, sigmaE = 0), nsim = 2)
  ff <- c(fbar(om$stock)[, "2025"]) / c(om$stock@refpts["Fmsy"])
  expect_equal(ff, rep(2, 2), tolerance = 0.05)
})

test_that("a selectivity change is used by the OM and the sampling", {
  skip_on_cran()
  dome <- selectivity("dnormal", Lpeak = 40, sd_left = 5, sd_right = 6)
  tr <- treatment("gear switch", selectivity = list(fleet = "linefish", from = 2015, sel = dome))
  om <- build_om(load_species("lutjanus_bohar"), reef_fishery(), nsim = 2, treatment = tr)
  expect_length(fleet_periods(om, "linefish"), 2)
  lfd <- sample_lfd(om, reef_oem(), 1)
  big <- function(y) {
    z <- lfd[lfd$source == "linefish" & lfd$year == y, ]
    sum(z$n[z$len_mm >= 600]) / sum(z$n)
  }
  expect_lt(big(2020), big(2012))     # dome gear: few large fish after the switch
})

test_that("length indicators respond to F and carry the truth", {
  skip_on_cran()
  sp <- one_gear_species(load_species("lutjanus_bohar"))
  oem <- oem_settings(one_gear_fishery(1), tail = 10, lfd_ess_age = 50, lfd_ess_len = 500)
  rt <- run_treatment(sp, one_gear_fishery(0.5),
                      treatment("hard", effort = list(fleet = "handline", from = 1969, mult = 4)),
                      nsim = 3, oem = oem, iters = 1:3)
  li <- length_indicators(rt, "Lutjanus bohar", iter = 1:3)
  m <- stats::aggregate(cbind(Lmean, Pmega, Z_lbar, Z_true) ~ run, li, mean)
  expect_lt(m$Lmean[m$run == "hard"], m$Lmean[m$run == "baseline"])
  expect_lt(m$Pmega[m$run == "hard"], m$Pmega[m$run == "baseline"])
  expect_gt(m$Z_lbar[m$run == "hard"], m$Z_lbar[m$run == "baseline"])
  expect_true(all(is.finite(li$Z_true)))
  expect_s3_class(plot_length_indicators(li), "ggplot")
  expect_s3_class(plot_indicators(rt$data$baseline, "Lutjanus bohar"), "ggplot")
})

test_that("planted artefacts are applied and listed in the answer key", {
  skip_on_cran()
  fy <- reef_fishery(creep = c(traps = 0, gillnet = 0, linefish = 0.03))
  oem <- reef_oem(fishery = fy, bio_age_bias = -0.2,
                  catch_bias = c(traps = 0, gillnet = 0, linefish = 0.3),
                  catch_bias_from = c(traps = 1971, gillnet = 1971, linefish = 2015),
                  cpue_q_step = list(gillnet = c(2020, 2)),
                  survey = survey_settings(q_step = c(2012, 0.5)))
  a <- planted_artefacts(oem, fy)
  expect_true(all(c("effort creep in CPUE", "catch bias", "CPUE catchability step",
                    "survey catchability step", "ageing bias") %in% a$artefact))
  om <- build_om(load_species("lutjanus_bohar"), fy, nsim = 2)
  d <- sample_data(om, oem, iters = 1)
  cc <- d$catch[d$catch$fleet == "linefish", ]
  expect_equal(mean(cc$catch_t[cc$year >= 2015] / cc$catch_true[cc$year >= 2015]),
               0.7, tolerance = 0.08)
  b <- d$bio[!is.na(d$bio$age) & d$bio$age_true >= 5, ]
  expect_lt(mean(b$age / b$age_true), 0.9)
  dir <- file.path(tempdir(), "floss_artefact_test")
  export_data(d, dir)
  expect_true(file.exists(file.path(dir, "answer_key", "planted_artefacts.csv")))
  expect_false(any(grepl("artefact", readLines(file.path(dir, "README.txt")))))
})

test_that("dynamics plots build", {
  skip_on_cran()
  sp <- load_species("siganus_sutor")
  om <- build_om(sp, reef_fishery(), nsim = 2)
  expect_s3_class(plot_eq(om), "ggplot")
  expect_s3_class(plot_srr(om, iters = 1:2), "ggplot")
  expect_s3_class(plot_ssb_vb(om), "ggplot")
  expect_s3_class(plot_eq_selectivity(sp, list(
    a = selectivity("logistic", L50 = 15, L95 = 20),
    b = selectivity("logistic", L50 = 25, L95 = 30))), "ggplot")
  expect_type(floss_tip("cpue_rho"), "character")
})
