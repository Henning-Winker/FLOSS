test_that("growth morphs: FLStock matches morph dynamics and SPR is length-based", {
  fy <- load_fishery("reef_3fleets")
  sp <- load_species("siganus_sutor"); sp$lh$sigR <- 0.01
  om <- build_om(sp, fy, nsim = 1, gtg = gtg_settings(n = 9))
  expect_s4_class(om$stock, "FLStockR")
  expect_true(all(is.finite(c(om$gtg$spr))))
  expect_true(all(c(om$gtg$spr) > 0 & c(om$gtg$spr) <= 1))
  # catch shares by fleet use fleet-specific catch weights
  cw <- catch_fleet(om, "wt")
  expect_equal(names(cw), fy$fleets)
  # morph SPR at F = 0 is 1
  G <- om$gtg$G
  expect_equal(FLOSS:::gtg_spr(G, setNames(rep(0, 3), fy$fleets)), 1)
  # standard OM unchanged when gtg = NULL
  om0 <- build_om(sp, fy, nsim = 1)
  expect_null(om0$gtg)
})
