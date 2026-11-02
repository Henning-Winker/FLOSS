test_that("effort trend reproduces the fudc shape", {
  yrs <- 1969:2025
  tr <- effort_trend(e0 = 0.05, fhi = 2.5, flo = 0.9, breaks = c(0.6, 0.8))
  v <- trend_values(tr, yrs)
  expect_length(v, length(yrs))
  expect_equal(v[1], 0.05)
  expect_equal(max(v), 2.5)
  expect_equal(tail(v, 1), 0.9)
})

test_that("effort is reproducible, shared and creep only affects efficiency", {
  fy <- reef_fishery(creep = 0.03)
  e1 <- effort_sim(fy, nsim = 4)
  e2 <- effort_sim(fy, nsim = 4)
  expect_identical(e1, e2)
  expect_equal(dims(e1$effort$traps)$iter, 4)
  e0 <- effort_sim(reef_fishery(), nsim = 4)
  expect_identical(c(e0$effort$traps), c(e1$effort$traps))
  expect_equal(c(e1$efficiency$traps[, "2025", , , , 1]), 1.03^(2025 - 1969))
  # caller's RNG state is untouched
  set.seed(1); a <- runif(1)
  set.seed(1); effort_sim(fy, nsim = 2); b <- runif(1)
  expect_identical(a, b)
})

test_that("build_om calibrates to MSY and partitions catch by fleet", {
  skip_on_cran()
  fy <- reef_fishery()
  sp <- load_species("siganus_sutor")
  om <- build_om(sp, fy, nsim = 5)
  expect_s3_class(om, "floss_om")
  expect_equal(an(FLBRP::refpts(om$brp)["msy", "yield"]), sp$msy,
               tolerance = 1e-3)
  expect_equal(dims(om$stock)$iter, 5)
  expect_equal(dims(om$stock)$minyear, 1971)
  cf <- catch_fleet(om)
  expect_named(cf, fy$fleets)
  expect_equal(c(Reduce(`+`, cf)), c(catch(om$stock)), tolerance = 1e-6)
  expect_true(all(c("Fmsy", "Bmsy", "MSY") %in% rownames(om$stock@refpts)))
})

test_that("build_oms shares effort and reuses existing OMs", {
  skip_on_cran()
  fy <- reef_fishery()
  oms <- build_oms(c("siganus_sutor", "lutjanus_bohar"), fy, nsim = 3)
  expect_identical(c(oms[[1]]$effort$linefish), c(oms[[2]]$effort$linefish))
  t0 <- oms$siganus_sutor$settings$built
  oms2 <- build_oms(c("siganus_sutor", "lethrinus_mahsena"), fy, nsim = 3,
                    oms = oms)
  expect_identical(oms2$siganus_sutor$settings$built, t0)
  expect_named(oms2, c("siganus_sutor", "lethrinus_mahsena"))
})
