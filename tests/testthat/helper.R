# Helper to load project functions for tests in a non-package project
load_project_R <- function() {
  # Locate project R/ directory robustly from tests/testthat
  r_dir <- NULL
  # Try testthat-relative path
  if (exists("test_path", where = asNamespace("testthat"), inherits = FALSE)) {
    r_dir <- testthat::test_path("..", "..", "R")
  }
  # Fallback to relative-from-cwd
  if (is.null(r_dir) || !dir.exists(r_dir)) {
    alt <- file.path("..", "..", "R")
    if (dir.exists(alt)) r_dir <- alt else r_dir <- "R"
  }
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) sys.source(f, envir = parent.frame())
}

load_project_R()
