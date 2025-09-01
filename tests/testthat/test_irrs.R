context("IRR computation")

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

test_that("IRR equals exp(coef) and has finite CIs", {
  load_project_R()
  df <- data.frame(
    site = paste0("S", 1:6), species = "sp", y = c(0,1,1,2,2,3),
    trap_nights = 10,
    type = factor(c("dry","wet","dry","wet","dry","wet"), levels = c("dry","wet")),
    leaf_litter_z = 0, cwd_z = 0, can_foliage_z = 0
  )
  m <- fit_species_glm(df, species = "sp")
  ti <- tidy_irrs(m)
  # pick a non-intercept term if present; otherwise intercept
  row <- ti[ti$term != "(Intercept)", , drop = FALSE]
  if (nrow(row) == 0) row <- ti[1, , drop = FALSE]
  coef_est <- coef(m$model)[row$term]
  expect_equal(as.numeric(row$irr), as.numeric(exp(coef_est)))
  expect_true(is.finite(row$lwr) && is.finite(row$upr))
})

