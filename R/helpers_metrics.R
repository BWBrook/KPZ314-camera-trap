# R/helpers_metrics.R

import::here(group_by, summarise, n_distinct, ungroup, tibble, .from = "dplyr")
import::here(pivot_wider, replace_na, .from = "tidyr")
import::here(diversity, estimateR, .from = "vegan")

#' Alpha‑diversity metrics per site
#'
#' @param df Data frame containing at minimum:
#'   * camera_site  (factor or character)
#'   * class_name   (factor or character)
#'   * event        (integer – independent events)
#'
#' @return tibble with one row per camera_site and columns:
#'   richness, shannon, simpson, pielou, chao1
#'
#' @export
calc_alpha <- function(df) {
  ## collapse to a site × species abundance table (event counts)
  mat <- df |>
    group_by(camera_site, class_name, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(camera_site, class_name) |>
    summarise(events = n(), .groups = "drop") |>
    pivot_wider(names_from = class_name,
                values_from = events,
                values_fill = 0L)

  site_id <- mat$camera_site
  mat_num <- as.matrix(mat[ , -1, drop = FALSE ])

  richness <- rowSums(mat_num > 0)
  shannon  <- diversity(mat_num, index = "shannon")
  simpson  <- diversity(mat_num, index = "simpson")
  pielou   <- shannon / log(pmax(richness, 1))
  chao1    <- unname(apply(mat_num, 1, function(v) estimateR(v)["S.chao1"]))

  tibble(
    camera_site = site_id,
    richness    = richness,
    shannon     = shannon,
    simpson     = simpson,
    pielou      = pielou,
    chao1       = chao1
  )
}

#' Community matrix (sites × species, abundance = event count)
#'
#' @param df Data frame with camera_site, class_name, event
#' @return matrix; row names = camera_site, col names = class_name
#' @export
build_comm_matrix <- function(df) {
  wide <- df |>
    group_by(camera_site, class_name, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(camera_site, class_name) |>
    summarise(events = n(), .groups = "drop") |>
    pivot_wider(
      names_from  = class_name,
      values_from = events,
      values_fill = 0L
    )

  mat <- as.matrix(wide[ , -1, drop = FALSE ])
  rownames(mat) <- wide$camera_site
  mat
}
