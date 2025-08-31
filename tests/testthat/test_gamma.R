context("gamma hill numbers")

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

test_that("gamma by habitat returns expected structure and coverage bounds", {
  load_project_R()
  # Two sites per habitat, simple events
  df <- data.frame(
    site = c(rep("D1",3), rep("D2",2), rep("W1",2), rep("W2",1)),
    common = c("A","A","B", "A","C", "A","B", "B"),
    event = c(1L,2L,1L, 1L,1L, 1L,1L, 1L),
    type  = c(rep("dry",5), rep("wet",3))
  )
  gh <- gamma_hill(df, group = "type", n_boot = 100L, seed = 1L)
  expect_true(all(gh$coverage >= 0 & gh$coverage <= 1, na.rm = TRUE))
  expect_true(all(gh$group == "type"))
  expect_true(all(c("dry","wet") %in% unique(gh$level)))
})

test_that("gamma all returns overall level", {
  load_project_R()
  df <- data.frame(
    site = c(rep("S1",3), rep("S2",1)),
    common = c("A","A","B","C"),
    event = c(1L,2L,1L,1L)
  )
  ga <- gamma_hill(df, group = "all", n_boot = 50L, seed = 1L)
  expect_true(all(ga$level == "overall"))
  expect_true(all(ga$q %in% c(0,1,2)))
})

