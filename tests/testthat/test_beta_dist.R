context("beta distances (Bray–Curtis and Jaccard)")

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

test_that("compute_dist returns symmetric dist objects with correct labels", {
  load_project_R()
  comm <- rbind(
    S1 = c(2,0,0),
    S2 = c(2,0,0),
    S3 = c(0,2,0),
    S4 = c(0,0,2)
  )
  colnames(comm) <- c("A","B","C")

  db <- compute_dist(comm, method = "bray")
  dj <- compute_dist(comm, method = "jaccard")

  expect_true(inherits(db, "dist"))
  expect_true(inherits(dj, "dist"))
  labs <- attr(db, "Labels"); if (is.null(labs)) labs <- attr(db, "labels")
  expect_true(all(labs %in% rownames(comm)))

  mb <- as.matrix(db); mj <- as.matrix(dj)
  expect_true(all(diag(mb) == 0))
  expect_true(all(diag(mj) == 0))
  expect_true(all(abs(mb - t(mb)) < 1e-12))
  expect_true(all(abs(mj - t(mj)) < 1e-12))
  # identical rows yield zero distance
  expect_equal(mb["S1","S2"], 0)
  expect_equal(mj["S1","S2"], 0)
})

