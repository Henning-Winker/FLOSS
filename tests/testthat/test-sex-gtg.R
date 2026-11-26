# Sex-structured OM with growth-type groups (gtg morphs) and morph-equilibrium reference points

lhF <- list(linf = 60, k = 0.15, t0 = -0.3, a = 1e-5, b = 3, L50 = 35, M = 0.3)
sel <- list(Trawl = selectivity("dnormal", Lpeak = 34, sd_left = 6, sd_right = 20))
g <- gtg_settings(n = 15, cv = 0.1, nsub = 12)

test_that("projection at Fmsy reaches morph Bmsy and MSY", {
  stk <- stock_sex(lh_sex(lhF, linf = 45, k = 0.25, L50 = 28), sel, tmax = 20, years = 2000:2002)
  Gs <- gtg_sex(stk, g)
  sr <- sr_sex(stk, s = 0.7, R0 = 1000, gtg = Gs)
  rp <- refpts_sex_gtg(Gs, sr)
  expect_true(all(c("Fmsy", "Bmsy", "MSY", "B0", "Fspr40") %in% dimnames(rp)$params))
  yrs <- 2000:2060
  om <- om_sex(stk, sr, fapic = rep(an(rp["Fmsy"]), length(yrs)), years = yrs, sigR = 0, gtg = Gs)
  expect_equal(sum(ssb(om$stock)[, "2060"]), an(rp["Bmsy"]), tolerance = 1e-3)
  expect_equal(sum(catch(om$stock)[, "2060"]), an(rp["MSY"]), tolerance = 1e-3)
})

test_that("identical sexes give half the one-sex morph spawning biomass per recruit", {
  stk1 <- stock_sex(lh_sex(lhF), sel, tmax = 20, years = 2000:2002)
  G1 <- gtg_sex(stk1, g)
  e <- eq_sex_gtg(G1, F = 0.3); e1 <- FLOSS:::gtg_eq(G1$F, c(Trawl = 0.3))
  expect_equal(e$sbpr, 0.5 * e1$phi, tolerance = 1e-12)
  expect_equal(e$ypr, e1$ypr, tolerance = 1e-12)
})

test_that("recruitment, truth and sexed length samples from the morphs", {
  stk <- stock_sex(lh_sex(lhF, linf = 45, k = 0.25, L50 = 28), sel, tmax = 20, years = 2000:2002)
  Gs <- gtg_sex(stk, g)
  sr <- sr_sex(stk, s = 0.7, R0 = 1000, gtg = Gs)
  fa <- fapic_fudc(2000:2020, nsim = 2, fref = refpts_sex_gtg(Gs, sr)["Fmsy"], seed = 3)
  om <- om_sex(stk, sr, fapic = fa, years = 2000:2020, sigR = 0.5, gtg = Gs)
  r <- stock.n(om$stock)[1]
  expect_equal(c(r[, , "F"]), c(r[, , "M"]))
  expect_true(all(c("F_Fmsy", "ssb_Bmsy", "spr") %in% names(om$truth)))
  lf <- lfd_sex(om, fleet = "Trawl", n = 200, sexed = TRUE, years = 2018:2020)
  expect_true(all(apply(lf, 2, sum) == 200))
  ls <- lfd_sex(om, survey_sel = selectivity("logistic", L50 = 12, L95 = 16), n = Inf, sexed = TRUE, years = 2020)
  expect_equal(sum(ls), 1, tolerance = 1e-8)
})
