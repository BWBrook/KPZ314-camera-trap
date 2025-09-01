context("GLM data preparation: scaling and baselines")

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

test_that("z-scored covariates approx mean 0, sd 1 and type baseline dry", {
  load_project_R()
  events <- data.frame(
    site = rep(paste0("S", 1:5), c(1,2,1,3,0)),
    common = c("sp1","sp1","sp1","sp1","sp1","sp1","sp1"),
    event = 1:7
  )
  site <- data.frame(
    site = paste0("S", 1:5),
    type = c("dry","wet","dry","wet","dry"),
    leaf_litter = c(1,2,3,4,5),
    cwd = c(0.1, 0.2, 0.2, 0.4, 0.5),
    can_foliage = c(10, 20, 10, 30, 30)
  )
  effort <- data.frame(site = paste0("S", 1:5), trap_nights = c(10,10,10,10,10))

  df <- prepare_species_glm_data(events, site, effort, species = "sp1",
                                 covars = c("type","leaf_litter","cwd","can_foliage"))
  expect_equal(levels(df$type), c("dry","wet"))
  expect_true(abs(mean(df$leaf_litter_z)) < 1e-8)
  expect_true(abs(sd(df$leaf_litter_z) - 1) < 1e-8)
})

