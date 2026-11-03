test_that("app files parse and run_app is available", {
  dir <- system.file("app", package = "FLOSS")
  expect_true(nzchar(dir))
  files <- c(file.path(dir, "app.R"), list.files(file.path(dir, "R"), full.names = TRUE))
  for (f in files) expect_silent(parse(f))
  expect_true(is.function(run_app))
})

test_that("app helpers adapt species to the fishery", {
  skip_if_not_installed("shiny")
  env <- new.env()
  sys.source(system.file("app", "R", "helpers.R", package = "FLOSS"), envir = env)
  sp <- load_species("lutjanus_bohar")
  one <- env$adapt_species(sp, one_gear_fishery(1), one_gear = TRUE)
  expect_equal(names(one$sel), "handline")
  reef <- env$adapt_species(sp, reef_fishery())
  expect_equal(names(reef$sel), c("traps", "gillnet", "linefish"))
})
