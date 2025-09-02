test_that("occupancy demo builds UF and fits when unmarked available", {
  load_project_R()
  testthat::skip_if_not_installed("unmarked")

  # simple 4-site, 5-day detection matrix
  hist <- matrix(
    c(1,0,0,0,0,
      0,1,0,0,0,
      0,0,1,0,0,
      0,0,0,1,1),
    nrow = 4, byrow = TRUE
  )
  rownames(hist) <- c("S1","S2","S3","S4")
  meta <- data.frame(site = rownames(hist), type = c("dry","dry","wet","wet"), stringsAsFactors = FALSE)

  umf <- occu_build_umf(hist, meta)
  fit <- fit_simple_occu(umf)
  td  <- fit$tidy
  expect_true(
    all(c("psi_dry","psi_wet","p") %in% td$param) ||
      all(c("psi_overall","p") %in% td$param)
  )
  expect_true(all(td$est > 0 & td$est < 1))
  expect_true(inherits(fit$fig, "ggplot"))
})
