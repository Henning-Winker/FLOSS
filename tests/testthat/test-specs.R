test_that("example library has the four reef species", {
  lib <- species_library()
  expect_true(all(c("siganus_sutor", "lethrinus_mahsena", "lutjanus_bohar",
                    "carcharhinus_melanopterus") %in% lib$id))
  sp <- load_species("lutjanus_bohar")
  expect_s3_class(sp, "floss_species")
  expect_equal(sp$msy, 600)
  expect_setequal(names(sp$sel), c("traps", "gillnet", "linefish"))
})

test_that("species and fishery survive a YAML round trip", {
  sp <- load_species("siganus_sutor")
  f <- tempfile(fileext = ".yaml")
  write_spec(sp, f)
  expect_identical(unclass_deep(read_spec(f)), unclass_deep(sp))

  fy <- reef_fishery(creep = c(traps = 0, gillnet = 0.01, linefish = 0.02))
  write_spec(fy, f)
  expect_identical(unclass_deep(read_spec(f)), unclass_deep(fy))
})

test_that("user species are saved, override examples and can be deleted", {
  withr_dir <- tempfile()
  old <- options(FLOSS.user_dir = withr_dir)
  on.exit(options(old))

  sp <- load_species("siganus_sutor")
  sp$msy <- 999
  save_species(sp)
  expect_error(save_species(sp), "exists")
  lib <- species_library()
  expect_equal(lib$source[lib$id == "siganus_sutor"], "user")
  expect_equal(load_species("siganus_sutor")$msy, 999)

  delete_species("siganus_sutor")
  expect_equal(load_species("siganus_sutor")$msy, 750)
  expect_error(delete_species("lutjanus_bohar"), "read-only")
})

test_that("invalid specs are rejected", {
  expect_error(selectivity("logistic", L50 = 20), "L95")
  sp <- load_species("siganus_sutor")
  lh <- sp$lh
  lh$s <- 1.5
  expect_error(species("x", lh = lh, msy = 1, sel = sp$sel,
                       relq = unlist(sp$relq)), "steepness")
  fy <- fishery(fleets = c("a", "b"), years = 2000:2020,
                trends = list(a = effort_trend(), b = effort_trend()))
  expect_error(validate_species(sp, fy), "no selectivity")
  expect_error(fishery(fleets = "a", years = 2000:2020,
                       trends = list(a = effort_trend(breaks = c(0.8, 0.2)))),
               "breaks")
})
