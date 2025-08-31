# R/helpers_effort.R — effort manifest and RAIs

import::here(select, mutate, transmute, rename, .from = "dplyr")
import::here(tibble, .from = "tibble")
import::here(cli_warn, .from = "cli")
import::here(left_join, group_by, summarise, n, ungroup, .from = "dplyr")

to_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  if (is.character(x)) {
    # Expect d/m/Y H:M:S
    return(as.Date(as.POSIXct(x, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")))
  }
  as.Date(NA)
}

clamp01 <- function(x) pmin(pmax(x, 0), 1)

#' Summarise sampling effort per site
#'
#' @param site_df Site metadata table, with columns site, op_days,
#'        first_image, last_image, days_with_event, blank, person, vehicle, total_images
#' @return tibble(site, trap_nights, days_with_event, blanks, person, vehicle, total_images, uptime_prop)
#' @export
summarise_effort <- function(site_df) {
  if (!all(c("site", "op_days") %in% names(site_df))) {
    rlang::abort("summarise_effort(): missing required columns 'site' and 'op_days'.")
  }

  first_d <- to_date_safe(site_df$first_image)
  last_d  <- to_date_safe(site_df$last_image)
  dur     <- as.numeric(last_d - first_d + 1)
  trap_n  <- as.numeric(site_df$op_days)
  uptime  <- ifelse(is.finite(dur) & is.finite(trap_n) & trap_n > 0,
                    clamp01(dur / trap_n),
                    NA_real_)

  out <- tibble::tibble(
    site           = site_df$site,
    trap_nights    = trap_n,
    days_with_event = site_df$days_with_event %||% NA_real_,
    blanks         = site_df$blank %||% NA_real_,
    person         = site_df$person %||% NA_real_,
    vehicle        = site_df$vehicle %||% NA_real_,
    total_images   = site_df$total_images %||% NA_real_,
    uptime_prop    = uptime
  )
  out
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Compute RAIs (per 100 trap-nights) by site and species
#'
#' @param events_df Data frame with site, common, event (independent events)
#' @param effort_df Output of summarise_effort()
#' @return tibble(site, common, events, trap_nights, rai_100)
#' @export
compute_rai <- function(events_df, effort_df) {
  if (!all(c("site", "common", "event") %in% names(events_df))) {
    rlang::abort("compute_rai(): events_df must contain 'site', 'common', 'event'.")
  }
  ev <- events_df |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site, common) |>
    summarise(events = n(), .groups = "drop")

  ef <- select(effort_df, site, trap_nights)
  out <- left_join(ev, ef, by = "site")

  zero_traps <- !is.na(out$trap_nights) & out$trap_nights == 0
  if (any(zero_traps)) cli_warn("compute_rai(): trap_nights = 0 for some sites; RAI set to NA.")

  out$rai_100 <- with(out, ifelse(is.na(trap_nights) | trap_nights == 0,
                                  NA_real_, 100 * events / pmax(trap_nights, 1)))
  out
}

