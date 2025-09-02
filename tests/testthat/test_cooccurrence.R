test_that("co-occurrence stats and z-matrix behave as expected", {
  load_project_R()
  # Build a tiny PA matrix with known margins
  # Sites S1..S6; species A,B,C
  pa <- rbind(
    c(1,1,0),
    c(1,0,1),
    c(1,1,1),
    c(0,1,0),
    c(0,0,1),
    c(1,0,0)
  )
  rownames(pa) <- paste0("S", 1:6)
  colnames(pa) <- c("A","B","C")

  st <- coocc_pair_stats(pa)
  expect_true(all(c("sp_i","sp_j","N","Xi","Xj","O","E","z","p_norm") %in% names(st)))
  # Check expectations for A,B
  N <- 6; Xi <- sum(pa[,"A"]); Xj <- sum(pa[,"B"]); O <- sum(pa[,"A"] & pa[,"B"])
  mu <- Xi * Xj / N
  var <- (Xi * Xj * (N - Xi) * (N - Xj)) / (N^2 * (N - 1))
  z  <- if (var > 0) (O - mu) / sqrt(var) else 0
  row_ab <- st[st$sp_i == "A" & st$sp_j == "B", , drop = FALSE]
  expect_equal(as.numeric(row_ab$E), mu)
  expect_equal(round(as.numeric(row_ab$z), 6), round(z, 6))

  # Symmetric z-matrix
  zm <- coocc_z_matrix(st)
  expect_true(is.matrix(zm))
  expect_true(all(dim(zm) == c(3,3)))
  expect_true(is.na(diag(zm)[1]))
  expect_equal(zm["A","B"], zm["B","A"])
})

