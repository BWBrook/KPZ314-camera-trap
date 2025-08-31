context("alpha hill numbers and bootstrap")

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

test_that("calc_alpha returns consistent Hill numbers", {
  load_project_R()
  # Construct events: site S1, species A with 3 events, B with 1
  df <- data.frame(
    site = c(rep("S1", 4)),
    common = c(rep("A", 3), "B"),
    event = c(1L, 2L, 3L, 1L)
  )
  out <- calc_alpha(df)
  expect_equal(nrow(out), 1)
  richness <- out$richness[1]
  shannon  <- out$shannon[1]
  invsimp  <- out$hill_q2[1]
  expect_equal(out$hill_q0[1], richness)
  expect_equal(out$hill_q1[1], exp(shannon))
  # inverse Simpson should match vegan's invsimpson
  expect_equal(invsimp, vegan::diversity(matrix(c(3,1), nrow = 1), index = "invsimpson"))
})

test_that("bootstrap_alpha CI contains point estimate", {
  load_project_R()
  df <- data.frame(
    site = rep("S1", 6),
    common = c("A","A","A","B","B","C"),
    event = c(1L,2L,3L,1L,2L,1L)
  )
  ci <- bootstrap_alpha(df, n_boot = 100L, seed = 1L)
  # ensure est lies within lwr..upr for all metrics
  expect_true(all(ci$lwr <= ci$est & ci$est <= ci$upr))
})

