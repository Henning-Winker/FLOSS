test_that("sex-structured equilibrium matches FLBRP for identical sexes and FLasher projections", {
  lhF <- list(linf = 60, k = 0.15, t0 = -0.3, a = 1e-5, b = 3, L50 = 35, M = 0.3)
  sel <- list(Trawl = selectivity("dnormal", Lpeak = 34, sd_left = 6, sd_right = 20))
  # identical sexes = one-sex FLBRP with half the recruitment
  stk1 <- stock_sex(lh_sex(lhF), sel, tmax = 20, years = 2000:2002)
  sr1 <- sr_sex(stk1, s = 0.7, R0 = 1000)
  rp1 <- refpts_sex(stk1, sr1)
  one <- stk1[, , "F"]
  b1 <- FLBRP::brp(FLBRP::FLBRP(one, sr = list(model = bevholt()$model,
    params = FLPar(a = an(params(sr1)["a"]) / 2, b = an(params(sr1)["b"])))))
  fb <- mean(harvest(one)[ac(range(one)["minfbar"]:range(one)["maxfbar"]), 1])
  expect_equal(an(rp1["Fmsy"]), an(FLBRP::refpts(b1)["msy", "harvest"]) / fb, tolerance = 1e-4)
  expect_equal(an(rp1["MSY"]), 2 * an(FLBRP::refpts(b1)["msy", "yield"]), tolerance = 1e-4)
  # dimorphic stock: FLasher holds the equilibrium at constant F
  lh <- lh_sex(lhF, linf = 45, k = 0.25, L50 = 28)
  stk <- stock_sex(lh, sel, tmax = 20, years = 2000:2002)
  sr <- sr_sex(stk, s = 0.7, R0 = 1000)
  st <- FLOSS:::init_stock_sex(stk, sr, F = 0.2, years = 2000:2030)
  o <- FLasher::fwd(st, sr = sr, control = FLasher::fwdControl(year = 2001:2030, quant = "effort", value = 0.2))
  e <- eq_sex(stk, sr, F = 0.2)
  expect_equal(sum(ssb(o)[, "2030"]), e$ssb, tolerance = 1e-6)
  expect_equal(sum(catch.n(o)[, "2030"] * catch.wt(o)[, "2030"]), e$yield, tolerance = 1e-6)
  # operating model and sexed length samples
  om <- om_sex(stk, sr, fapic = rep(0.2, 3), sigR = 0)
  l <- lfd_sex(om, survey_sel = selectivity("logistic", L50 = 12, L95 = 16), n = 100, sexed = TRUE)
  expect_equal(dimnames(l)$unit, c("F", "M"))
  expect_true(all(apply(l, 2, sum) == 100))
})
