test_that("effort in fishing days, aggregation, raising and FLR export", {
  skip_on_cran()
  fy <- load_fishery()
  oms <- build_oms(c("lutjanus_bohar", "siganus_sutor"), fy, nsim = 2)
  d <- sample_data(oms, reef_oem(), iters = 1)
  expect_true(all(d$effort$effort == round(d$effort$effort)))
  expect_gt(max(d$effort$effort), 100)
  ef <- effort_table(d, "aggregate")
  expect_equal(sum(ef$effort), sum(d$effort$effort))
  ct <- catch_table(d, by = "combined")
  expect_equal(sum(ct$catch_t), sum(d$catch$catch_t))
  nc <- nominal_cpue(d, "Lutjanus bohar", "aggregate")
  expect_equal(unique(nc$source), "all fleets")
  expect_equal(mean(nc$index), 1)
  # raised length samples weigh as much as the catch
  lr <- lfd_aggregate(d, "Lutjanus bohar", "fleet", "catch")
  lh <- d$lhpars[d$lhpars$species == "Lutjanus bohar", ]
  bw <- FLOSS:::lfd_bin_mm(d, "Lutjanus bohar") / 10
  z <- lr[lr$year == 2025 & lr$source == "linefish", ]
  w <- sum(z$n * lh$a * (z$len_mm / 10 + bw / 2)^lh$b) / 1e6
  cc <- d$catch$catch_t[d$catch$species == "Lutjanus bohar" & d$catch$year == 2025 & d$catch$fleet == "linefish"]
  expect_equal(w, cc, tolerance = 0.01)
  la <- lfd_aggregate(d, "Lutjanus bohar", "aggregate", "sample")
  expect_equal(sum(la$n), sum(d$lfd$n[d$lfd$species == "Lutjanus bohar" & d$lfd$source != "survey"]))
  expect_s3_class(plot_aggregate(d, "Lutjanus bohar", what = "cpue", by = "combined"), "ggplot")
  expect_s3_class(plot_bio(d, "Lutjanus bohar", by = "season"), "ggplot")
  f <- as_flr(d, "Lutjanus bohar", 1)
  expect_s4_class(f$catch, "FLQuants")
  expect_s4_class(f$indices, "FLIndices")
  expect_true(all(is.finite(c(index.var(f$indices[["linefish"]])))))
  expect_equal(dimnames(f$lfd[[1]])[[1]][1:2] |> as.numeric() |> diff(), bw)
  expect_equal(as.numeric(f$lhpar["a"]), lh$a / 1000)
  expect_equal(sum(f$catch_by_gear), 1)
  tmp <- tempfile()
  export_data(d, tmp, format = c("csv", "rds"), aggregate = TRUE)
  expect_true(all(c("floss_flr.rds", "cpue_nominal.csv", "lfd_aggregate.csv") %in%
                    list.files(file.path(tmp, "iter_01"))))
  expect_named(readRDS(file.path(tmp, "iter_01", "floss_flr.rds")), c("Lutjanus bohar", "Siganus sutor"))
})
