# helpers_extend.R — additional KPI functions for KPZ314 prac
# -----------------------------------------------------------------------------
# These lightweight helpers slot into _targets.R for the “stretch” analyses:
#   * calc_rarefaction():  ggplot of rarefaction curves for a community matrix
#   * permanova_region():  PERMANOVA on Bray–Curtis vs a `region` factor
#   * species_trend():     Event‑count table for focal species across sites
# 
# All follow the explicit‑import rule (import::from), no side effects.
# -----------------------------------------------------------------------------

import::from("vegan", rarefy, adonis2)
import::from("ggplot2", ggplot, aes, geom_line, theme_minimal, labs)
import::from("dplyr", bind_rows, distinct, filter, count)
import::from("tidyr", pivot_wider)

#' Rarefaction curves (one per camera site, rows = sites)
#'
#' @param comm_mat matrix of event counts, rows = sites, cols = species
#' @param step     positive integer increment for subsampling
#' @return         ggplot object
#' @export
calc_rarefaction <- function(comm_mat, step = 1) {
  ## 1. ensure rows = sites, integer counts
  comm_mat <- as.matrix(comm_mat)
  storage.mode(comm_mat) <- "integer"

  curves <- vector("list", nrow(comm_mat))

  for (i in seq_len(nrow(comm_mat))) {
    counts <- comm_mat[i, ]
    n_tot  <- sum(counts)

    # skip sites with fewer than `step` events
    if (n_tot < step) next

    subs  <- seq(step, n_tot, by = step)
    s_hat <- rarefy(counts, sample = subs)

    curves[[i]] <- data.frame(
      site  = rownames(comm_mat)[i],
      reads = subs,
      S     = as.numeric(s_hat)
    )
  }

  df <- bind_rows(curves)

  ggplot(df, aes(reads, S, group = site)) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Rarefaction curves",
      x     = "Number of events (subsample)",
      y     = expression("Expected species richness ("*S[hat]*")")
    )
}

#' PERMANOVA for wet vs dry (or any) regions
#'
#' @param dist_obj  vegan::vegdist object (Bray–Curtis)
#' @param meta_df   data frame with columns camera_site and region
#' @return          adonis2 result (data.frame)
#' @export
permanova_region <- function(dist_obj, meta_df) {
  # adonis2 expects one row per site in the metadata
  meta_unique <- distinct(meta_df, camera_site, region)
  adonis2(dist_obj ~ region, data = meta_unique)
}

#' Event counts for focal species across sites
#'
#' @param df           joined data frame with columns camera_site, class_name, event
#' @param species_vec  character vector of focal species (class_name values)
#' @return             wide table: camera_site × species, values = event counts
#' @export
species_trend <- function(df, species_vec) {
  df |>
    filter(class_name %in% species_vec) |>
    count(camera_site, class_name, name = "events") |>
    pivot_wider(names_from = class_name, values_from = events, values_fill = 0)
}
