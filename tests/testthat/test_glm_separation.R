context("GLM separation guard for type")

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

test_that("type dropped when all events occur in one level", {
  load_project_R()
  df <- data.frame(
    site = paste0("S", 1:6), species = "sp",
    y = c(2,1,0,0,0,0),
    trap_nights = 10,
    type = factor(c("dry","dry","wet","wet","wet","wet"), levels = c("dry","wet")),
    leaf_litter_z = 0, cwd_z = 0, can_foliage_z = 0
  )
  m <- fit_species_glm(df, species = "sp")
  dg <- diagnostics_table(m)
  expect_true(grepl("type dropped", dg$note))
})

