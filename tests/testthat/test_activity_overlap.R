test_that("activity overlap responds to separation and returns expected structure", {
  load_project_R()
  # Synthetic events: species A mostly day, B mostly night
  make_dt <- function(hours_vec) as.POSIXct(sprintf("2025-01-01 %02d:00:00", hours_vec), tz = "Australia/Hobart")
  df <- data.frame(
    site = rep(paste0("S", 1:6), each = 4),
    common = rep(c("A","B"), each = 12),
    event = rep(1:24),
    datetime = c(make_dt(sample(8:16, 12, replace = TRUE)),
                 make_dt(sample(c(0:5,20:23), 12, replace = TRUE))),
    type = rep(rep(c("dry","wet"), each = 12), 1),
    stringsAsFactors = FALSE
  )

  ov <- overlap_by_pair(df, c("A","B"), groups = c("dry","wet"), n_boot = 50L, seed = 1L)
  expect_true(all(c("species1","species2","type","n1","n2","estimator","delta_hat","lwr","upr") %in% names(ov)))
  # plot returns ggplot
  p <- plot_activity_overlap(df, c("A","B"), groups = c("dry","wet"))
  expect_true(inherits(p, "ggplot"))
})

