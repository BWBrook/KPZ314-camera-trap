context("beta inference: NMDS + PERMANOVA + dispersion")

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

test_that("permanova_and_dispersion returns expected columns and ranges", {
  load_project_R()
  # construct a 6-site matrix with 3 per group
  comm <- rbind(
    D1 = c(5,0,0,0),
    D2 = c(4,1,0,0),
    D3 = c(6,0,0,0),
    W1 = c(0,0,5,0),
    W2 = c(0,0,4,1),
    W3 = c(0,0,6,0)
  )
  colnames(comm) <- c("A","B","C","D")
  site <- data.frame(site = rownames(comm), type = rep(c("dry","wet"), each = 3))

  db <- compute_dist(comm, method = "bray")
  dj <- compute_dist(comm, method = "jaccard")

  tb <- permanova_and_dispersion(db, site, group_col = "type", n_perm = 199L, seed = 1L)
  tj <- permanova_and_dispersion(dj, site, group_col = "type", n_perm = 199L, seed = 1L)

  expect_true(all(c("distance","term","df","pseudo_F","R2","p_perm","disp_F","disp_p","note") %in% names(tb)))
  expect_true(all(tb$R2 >= 0 & tb$R2 <= 1))
  expect_true(all(tb$p_perm >= 0 & tb$p_perm <= 1))
  expect_true(all(tj$R2 >= 0 & tj$R2 <= 1))
  expect_true(all(tj$p_perm >= 0 & tj$p_perm <= 1))
})

test_that("nmds_with_hulls returns ggplot and stats with stress", {
  load_project_R()
  comm <- rbind(
    D1 = c(5,0,0,0),
    D2 = c(4,1,0,0),
    D3 = c(6,0,0,0),
    W1 = c(0,0,5,0),
    W2 = c(0,0,4,1),
    W3 = c(0,0,6,0)
  )
  colnames(comm) <- c("A","B","C","D")
  meta <- data.frame(site = rownames(comm), type = rep(c("dry","wet"), each = 3))
  db <- compute_dist(comm, method = "bray")
  nm <- suppressWarnings(nmds_with_hulls(db, meta, group_col = "type", k = 2L, seed = 1L))
  expect_true(inherits(nm$fig, "ggplot"))
  expect_true(is.data.frame(nm$stats))
  expect_true(is.numeric(nm$stats$stress))
  # has a polygon layer when >= 3 per group
  has_poly <- any(vapply(nm$fig$layers, function(ly) inherits(ly$geom, "GeomPolygon"), logical(1)))
  expect_true(has_poly)
})
