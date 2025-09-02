test_that("unknown_animal sensitivity behaves as expected", {
  load_project_R()

  # small synthetic events with an unknown_animal record
  events <- data.frame(
    site = c("A","A","B","B","B","C","C"),
    common = c("sp1","unknown_animal","sp1","sp2","unknown_animal","sp2","sp3"),
    event = c(1L,2L,1L,2L,3L,1L,2L),
    datetime = as.POSIXct("2025-01-01 00:00:00", tz = "UTC") + 0:6,
    stringsAsFactors = FALSE
  )

  a_all <- calc_alpha(events)
  a_no  <- calc_alpha(events_without_unknown(events))
  del   <- alpha_delta_unknown(a_all, a_no)
  # richness deltas must be <= 0
  rchg <- del[del$metric == "richness", , drop = FALSE]
  expect_true(all(rchg$delta <= 1e-12))

  c_all <- build_comm_matrix(events)
  c_no  <- build_comm_matrix(events_without_unknown(events))
  bd <- beta_delta_summary(c_all, c_no)
  expect_true(all(bd$distance %in% c("bray","jaccard")))
  expect_true(all(is.finite(bd$r_upper)))
  expect_true(all(bd$r_upper >= -1 & bd$r_upper <= 1))
})

