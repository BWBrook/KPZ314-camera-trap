context("rarefaction plots")

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

test_that("rarefaction functions return ggplot objects", {
  load_project_R()
  # site × species counts
  comm <- rbind(c(3,1,0), c(0,2,2))
  colnames(comm) <- c("A","B","C")
  rownames(comm) <- c("S1","S2")
  site <- data.frame(site = c("S1","S2"), type = c("dry","wet"))
  p1 <- rarefaction_site(comm, step = 1L)
  p2 <- rarefaction_by_group(comm, site, group_col = "type", step = 1L)
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
})

