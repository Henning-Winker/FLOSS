oms_small <- function() {
  if (is.null(the$oms)) the$oms <- build_oms(c("siganus_sutor", "lutjanus_bohar"),
                                             reef_fishery(), nsim = 3)
  the$oms
}
the <- new.env()

test_that("oem settings take per-fleet values and tails", {
  oem <- reef_oem()
  expect_s3_class(oem, "floss_oem")
  expect_equal(oem$years$traps, c(2021L, 2025L))
  expect_equal(oem$years$linefish, c(2011L, 2025L))
  expect_equal(oem$cpue$rho$linefish, 0.6)
  expect_error(oem_settings(reef_fishery(), tail = c(traps = 5)), "no value")
})

test_that("AR1 observation error has the requested autocorrelation", {
  set.seed(1)
  x <- FLQuant(1, dimnames = list(year = 1:2000))
  o <- ar1_obs(x, cv = 0.3, rho = 0.7)
  e <- log(c(o))
  expect_equal(stats::acf(e, plot = FALSE)$acf[2], 0.7, tolerance = 0.05)
  expect_equal(mean(c(o)), 1, tolerance = 0.03)
})

test_that("sample_data returns all streams within the sampling windows", {
  skip_on_cran()
  dat <- sample_data(oms_small(), reef_oem(), iters = 1:2, seed = 1)
  expect_s3_class(dat, "floss_data")
  expect_setequal(unique(dat$catch$iter), 1:2)
  expect_true(all(dat$cpue$year[dat$cpue$source == "traps"] >= 2021))
  expect_true(all(dat$lfd$year[dat$lfd$source == "linefish"] >= 2011))
  # two-stage LFD: ess_len fish per year and source
  n <- aggregate(n ~ species + iter + year + source, dat$lfd, sum)
  expect_true(all(n$n[n$source == "linefish"] == 500))
  expect_true(all(n$n[n$source == "traps"] == 300))
  # biological samples: one row per fish, ages partly missing
  b <- dat$bio
  expect_equal(nrow(b[b$species == "Siganus sutor" & b$iter == 1 &
                        b$source == "traps" & b$year == 2025, ]), 100)
  expect_true(all(b$mature %in% 0:1))
  expect_true(any(is.na(b$age)) && any(!is.na(b$age)))
  # reproducible
  expect_identical(sample_data(oms_small(), reef_oem(), iters = 1:2, seed = 1)$bio,
                   dat$bio)
})

test_that("stratified ageing caps otoliths per length bin", {
  skip_on_cran()
  oem <- reef_oem(bio_age_scheme = "stratified", bio_n_per_bin = 2)
  b <- sample_bio(oms_small()$siganus_sutor, oem, 1)
  b <- b[!is.na(b$age), ]
  per_bin <- table(b$year, b$source, floor(b$len_mm / 10))
  expect_lte(max(per_bin), 2)
})

test_that("ALK, catch-at-age and export work on observed data", {
  skip_on_cran()
  dat <- sample_data(oms_small(), reef_oem(), iters = 1, seed = 2)
  k <- alk_obs(dat, "Lutjanus bohar", prop = TRUE)
  expect_true(all(abs(rowSums(k, na.rm = TRUE)[!is.na(rowSums(k))] - 1) < 1e-9))
  caa <- caa_obs(dat, "Lutjanus bohar", "linefish")
  lf <- aggregate(n ~ year, dat$lfd[dat$lfd$species == "Lutjanus bohar" &
                                      dat$lfd$source == "linefish", ], sum)
  tot <- aggregate(n ~ year, caa, sum)
  un <- attr(caa, "unassigned")
  expect_equal(tot$n, lf$n * (1 - un[as.character(lf$year)]), tolerance = 1e-6,
               ignore_attr = TRUE)
  dir <- file.path(tempdir(), "floss_export_test")
  export_data(dat, dir, give = "data")
  expect_true(file.exists(file.path(dir, "iter_01", "bio.csv")))
  expect_false(file.exists(file.path(dir, "iter_01", "lhpars.csv")))
  expect_true(file.exists(file.path(dir, "answer_key", "truth.csv")))
  d1 <- read_dataset(file.path(dir, "iter_01"))
  expect_false("age_true" %in% names(d1$bio))
  expect_false("catch_true" %in% names(d1$catch))
})

test_that("true status summary and truth export", {
  skip_on_cran()
  dat <- sample_data(oms_small(), reef_oem(), iters = 1:2, seed = 3)
  st <- status_table(dat)
  expect_equal(nrow(st), 4)
  expect_true(all(st$status %in% c("Healthy", "Overfishing", "Overfished",
                                   "Overfished and overfishing")))
  om_st <- status_table(oms_small()$siganus_sutor, iters = 1:2)
  expect_equal(st$B_Bmsy[st$species == "Siganus sutor"], om_st$B_Bmsy)
  s <- status_table(dat, summary = TRUE)
  expect_equal(nrow(s), 2)
  dir <- file.path(tempdir(), "floss_truth_test")
  export_truth(dat, dir, format = c("csv", "xlsx"))
  expect_true(all(file.exists(file.path(dir, c("truth_timeseries.csv",
    "refpts.csv", "status.csv", "truth.xlsx")))))
  expect_s3_class(plot_truth(dat), "ggplot")
  expect_s3_class(plot_kobe(dat), "ggplot")
})
