test_that("detection history builds correctly and ψ̂(t) is non-decreasing", {
  load_project_R()

  # synthetic site metadata
  site_df <- data.frame(
    site = c("S1","S2","S3"),
    type = c("dry","wet","dry"),
    first_image = as.Date(c("2025-01-01","2025-01-01","2025-01-01")),
    op_days = c(10L, 10L, 10L),
    stringsAsFactors = FALSE
  )
  # synthetic events for species sp1
  events <- data.frame(
    site = c("S1","S1","S2","S3"),
    common = c("sp1","sp1","sp1","sp1"),
    event = c(1L,2L,1L,1L),
    datetime = as.POSIXct(c("2025-01-01 10:00:00","2025-01-03 09:00:00","2025-01-02 12:00:00","2025-01-05 08:00:00"), tz = "UTC"),
    stringsAsFactors = FALSE
  )

  dh <- detection_history_for_species(events, site_df, species = "sp1", max_days = 7L)
  expect_equal(dim(dh$mat), c(3,7))
  expect_equal(rownames(dh$mat), site_df$site)

  curves <- naive_psi_curves(dh$mat, dh$meta)
  # check non-decreasing within each type
  by_grp <- split(curves, interaction(curves$species, curves$type))
  for (df in by_grp) {
    diffs <- diff(df$psi_naive)
    expect_true(all(diffs >= -1e-12))
  }
})

