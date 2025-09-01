context("GLM family selection and diagnostics")

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

test_that("Poisson chosen under low dispersion; NB considered under high dispersion", {
  load_project_R()
  # Construct a small df with low dispersion (y ~ Poisson-like)
  df1 <- data.frame(
    site = paste0("S", 1:6), species = "sp", y = c(0,1,1,2,2,3),
    trap_nights = 10,
    type = factor(c("dry","wet","dry","wet","dry","wet"), levels = c("dry","wet")),
    leaf_litter_z = 0, cwd_z = 0, can_foliage_z = 0
  )
  m1 <- fit_species_glm(df1, species = "sp")
  expect_equal(m1$family, "poisson")

  # High overdispersion: introduce extreme variability
  df2 <- data.frame(
    site = paste0("S", 1:10), species = "sp", y = c(0,0,0,1,2,3,10,20,30,40),
    trap_nights = 10,
    type = factor(rep(c("dry","wet"), each = 5), levels = c("dry","wet")),
    leaf_litter_z = 0, cwd_z = 0, can_foliage_z = 0
  )
  m2 <- fit_species_glm(df2, species = "sp")
  expect_true(m2$family %in% c("negbin","quasipoisson"))
})

