# R/helpers_sensitivity.R — gap sensitivity diagnostics

import::here(group_by, summarise, n_distinct, ungroup, left_join, mutate, select, .from = "dplyr")
import::here(tibble, .from = "tibble")
import::here(pivot_longer, .from = "tidyr")
import::here(ggplot, aes, geom_line, facet_wrap, theme_minimal, labs, .from = "ggplot2")

#' Summarise sensitivity of metrics to the event gap
#'
#' @param events_list list of event-labelled data frames (each with attr 'gap_min')
#' @param effort_df effort manifest (summarise_effort)
#' @return list(metrics = tibble, fig = ggplot, deltas = tibble)
#' @export
summarise_gap_sensitivity <- function(events_list, effort_df) {
  if (!length(events_list)) {
    return(list(metrics = tibble(), fig = ggplot2::ggplot(), deltas = tibble()))
  }

  metrics_list <- lapply(events_list, function(ev) {
    gap <- attr(ev, "gap_min") %||% NA_integer_

    ev_unique <- ev |>
      group_by(site, common, event) |>
      summarise(n = 1L, .groups = "drop")

    by_site <- ev_unique |>
      group_by(site) |>
      summarise(
        richness = n_distinct(common),
        total_events = n_distinct(event),
        .groups = "drop"
      )

    # RAI per site (mean across species)
    rai_tbl <- compute_rai(ev, effort_df) |>
      group_by(site) |>
      summarise(mean_rai_100 = mean(rai_100, na.rm = TRUE), .groups = "drop")

    left_join(by_site, rai_tbl, by = "site") |>
      mutate(gap_min = as.integer(gap))
  })

  metrics <- do.call(rbind, metrics_list)

  # Relative deltas vs 5 min baseline
  baseline <- metrics |>
    group_by(site) |>
    summarise(
      richness_5 = richness[gap_min == 5][1],
      total_events_5 = total_events[gap_min == 5][1],
      mean_rai_100_5 = mean_rai_100[gap_min == 5][1],
      .groups = "drop"
    )

  deltas <- left_join(metrics, baseline, by = "site") |>
    mutate(
      richness = ifelse(is.finite(richness_5) & richness_5 > 0, richness / richness_5 - 1, NA_real_),
      total_events = ifelse(is.finite(total_events_5) & total_events_5 > 0, total_events / total_events_5 - 1, NA_real_),
      mean_rai_100 = ifelse(is.finite(mean_rai_100_5) & mean_rai_100_5 > 0, mean_rai_100 / mean_rai_100_5 - 1, NA_real_)
    ) |>
    select(site, gap_min, richness, total_events, mean_rai_100) |>
    pivot_longer(-c(site, gap_min), names_to = "metric", values_to = "rel_change_vs_5min")

  # Figure: metrics vs gap
  plot_df <- metrics |>
    pivot_longer(-c(site, gap_min), names_to = "metric", values_to = "value")

  fig <- ggplot(plot_df, aes(x = gap_min, y = value, group = site)) +
    geom_line(alpha = 0.4) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal() +
    labs(
      title = "Gap-sensitivity diagnostics",
      x = "Gap (minutes)", y = "Value"
    )

  list(metrics = metrics, fig = fig, deltas = deltas)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

