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

test_that("gtg_settings partitions the CV of Linf exactly", {
  s <- gtg_settings(n = 15, cv = 0.1)
  expect_equal(s$cv_res, 0.1 / sqrt(15))
  expect_equal(s$cv_between^2 + s$cv_res^2, 0.1^2)
  expect_error(gtg_settings(cv = 0.1, cv_res = 0.2))
  sp <- load_species("siganus_sutor")
  G <- FLOSS:::gtg_setup(sp, "traps", s)
  # weighted SD of the morph Linf grid equals the between-morph CV (the grid is rescaled for truncation)
  m <- sum(G$recP * G$linf_g)
  expect_equal(sqrt(sum(G$recP * (G$linf_g - m)^2)) / sp$lh$linf, s$cv_between)
})

test_that("true catch-at-length matches the morph keys and catch at age", {
  fy <- load_fishery("reef_3fleets")
  sp <- load_species("siganus_sutor"); sp$lh$sigR <- 0.01
  om <- build_om(sp, fy, nsim = 2, gtg = gtg_settings(n = 9))
  lq <- gtg_lfq(om, seasons = 4, iters = 2)
  expect_equal(names(lq), fy$fleets)
  expect_equal(dim(lq[[1]])[c(4, 6)], c(4L, 1L))
  # all of each fleet's catch numbers within the length range are allocated to lengths
  cn <- catch_fleet(om, "n")
  y <- dimnames(om$stock)$year[10]
  for (f in fy$fleets)
    expect_equal(sum(lq[[f]][, y]), sum(cn[[f]][, y, , , , 2]), tolerance = 0.02)
})
