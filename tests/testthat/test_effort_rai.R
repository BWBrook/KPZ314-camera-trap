context("effort and RAI")

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

test_that("effort summary computes uptime and mapping", {
  load_project_R()
  site <- data.frame(
    site = c("S1", "S2"),
    op_days = c(10, 0),
    first_image = c("01/01/2020 00:00:00", "01/01/2020 00:00:00"),
    last_image  = c("10/01/2020 23:59:59", "01/01/2020 00:00:00"),
    days_with_event = c(5, 0), blank = 0, person = 0, vehicle = 0, total_images = 0
  )
  ef <- summarise_effort(site)
  expect_equal(names(ef), c("site","trap_nights","days_with_event","blanks","person","vehicle","total_images","uptime_prop"))
  expect_equal(ef$trap_nights[1], 10)
  expect_true(is.na(ef$uptime_prop[2]))  # divide by zero -> NA
})

test_that("RAI per 100 trap-nights computed and zeros yield NA", {
  load_project_R()
  events <- data.frame(
    site = c("S1","S1","S2","S2"),
    common = c("a","a","b","b"),
    event = c(1L, 2L, 1L, 2L)
  )
  effort <- data.frame(site = c("S1","S2"), trap_nights = c(10, 0))
  expect_warning(
    out <- compute_rai(events, effort),
    regexp = "trap_nights = 0"
  )
  s1 <- subset(out, site == "S1")
  expect_equal(s1$events, 2)
  expect_equal(s1$rai_100, 100 * 2/10)
  s2 <- subset(out, site == "S2")
  expect_true(all(is.na(s2$rai_100)))
})
