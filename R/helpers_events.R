# R/helpers_events.R — event building and gap parameterisation

import::here(group_by, arrange, mutate, ungroup, select, .from = "dplyr")
import::here(tibble, .from = "tibble")
import::here(abort, syms, .from = "rlang")

#' Build independent events based on a time gap threshold
#'
#' @param df Data frame with at least columns: site, common, datetime
#' @param gap_min Integer minutes; detections within gap_min remain same event
#' @param group Character vector of grouping columns (default c("site","common"))
#'
#' @return Same rows with integer column `event` overwritten/added.
#'         Attribute `gap_min` is set on the returned data frame.
#'
#' @details Invariants:
#' * Within each group, event IDs start at 1 and only increment when the
#'   difference to the previous detection exceeds `gap_min` minutes.
#' * Rows are retained. The function requires `datetime` to be present and
#'   sorted within each group; otherwise it aborts.
#'
#' @export
build_events <- function(df, gap_min = 5L, group = c("site", "common")) {
  if (!all(c("datetime", group) %in% names(df))) {
    abort("build_events(): missing required columns 'datetime', 'site', or 'common'.")
  }
  if (!is.numeric(gap_min) || length(gap_min) != 1L || gap_min < 0) {
    abort("build_events(): 'gap_min' must be a single non-negative number.")
  }

  # Ensure POSIXct datetime; abort on NA
  dt <- df$datetime
  if (inherits(dt, "character")) {
    dt <- as.POSIXct(dt, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")
  }
  if (!inherits(dt, "POSIXct")) {
    abort("build_events(): 'datetime' must be POSIXct or character in d/m/Y H:M:S format.")
  }
  if (anyNA(dt)) abort("build_events(): 'datetime' contains NA; cannot build events.")

  df$datetime <- dt

  # Check sortedness within groups
  is_sorted <- function(x) {
    # non-decreasing
    all(!is.na(x)) && all(order(x) == seq_along(x))
  }

  # split by group to validate ordering
  grp_split <- split(seq_len(nrow(df)), f = interaction(df[group], drop = TRUE))
  bad <- vapply(grp_split, function(idx) !is_sorted(df$datetime[idx]), logical(1))
  if (any(bad)) {
    abort("build_events(): 'datetime' must be sorted within groups; found unsorted groups.")
  }

  # Compute events within groups
  out <- df |>
    group_by(!!!syms(group)) |>
    arrange(datetime, .by_group = TRUE) |>
    mutate(
      .delta_min = as.numeric(difftime(datetime, dplyr::lag(datetime), units = "mins")),
      event = {
        # first row in each group is event 1; increment when gap exceeded
        inc <- ifelse(is.na(.delta_min) | .delta_min <= gap_min, 0L, 1L)
        cumsum(replace(inc, 1L, 1L))
      }
    ) |>
    ungroup() |>
    select(-.delta_min)

  attr(out, "gap_min") <- as.integer(gap_min)
  out
}
