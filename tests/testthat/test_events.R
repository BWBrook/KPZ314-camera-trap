context("build_events")

# Ensure project helpers are loaded for non-package tests
load_project_R <- function() {
  r_dir <- NULL
  if (exists("test_path", where = asNamespace("testthat"), inherits = FALSE)) {
    r_dir <- testthat::test_path("..", "..", "R")
  }
  if (is.null(r_dir) || !dir.exists(r_dir)) {
    alt <- file.path("..", "..", "R")
    if (dir.exists(alt)) r_dir <- alt else r_dir <- "R"
  }
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) sys.source(f, envir = parent.frame())
}
load_project_R()

test_that("gap splitting works at 5 minutes", {
  load_project_R()
  # Gaps: 1, 4, 6 minutes -> expect last to split at 5
  df <- data.frame(
    site = "S1",
    common = "sp",
    datetime = as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + c(0, 60, 300, 660),
    stringsAsFactors = FALSE
  )
  ev <- build_events(df, gap_min = 5L)
  expect_equal(ev$event, c(1L, 1L, 1L, 2L))
})

test_that("unsorted times abort", {
  load_project_R()
  df <- data.frame(
    site = c("S1", "S1"),
    common = c("sp", "sp"),
    datetime = as.POSIXct(c("2020-01-01 00:01:00", "2020-01-01 00:00:00"), tz = "UTC")
  )
  expect_error(build_events(df, gap_min = 5L))
})

test_that("NA datetime aborts", {
  load_project_R()
  df <- data.frame(
    site = c("S1", "S1"),
    common = c("sp", "sp"),
    datetime = as.POSIXct(c(NA, "2020-01-01 00:00:00"), tz = "UTC")
  )
  expect_error(build_events(df, gap_min = 5L))
})

test_that("mixed sites and species handled", {
  load_project_R()
  df <- data.frame(
    site = c("S1", "S1", "S2", "S2"),
    common = c("a", "a", "b", "b"),
    datetime = as.POSIXct("2020-01-01 00:00:00", tz = "UTC") + c(0, 60, 0, 600)
  )
  ev <- build_events(df, gap_min = 5L)
  # Each group starts at 1 and increments appropriately
  expect_true(all(ev$event[c(1,2)] == c(1L,1L)))
  expect_true(all(ev$event[c(3,4)] == c(1L,2L)))
})
