# helpers_extend.R — additional KPI functions for KPZ314 prac
# -----------------------------------------------------------------------------
# These lightweight helpers slot into _targets.R for the “stretch” analyses:
#   * calc_rarefaction():  ggplot of rarefaction curves for a community matrix
#   * permanova_region():  PERMANOVA on Bray–Curtis vs a `type` factor
#   * species_trend():     Event‑count table for focal species across sites
# 
# All follow the explicit‑import rule (import::from), no side effects.
# -----------------------------------------------------------------------------

import::here(rarefy, adonis2, .from = "vegan")
import::here(ggplot, aes, geom_line, theme_minimal, labs, .from = "ggplot2")
import::here(bind_rows, distinct, filter, count, .from = "dplyr")
import::here(pivot_wider, .from = "tidyr")

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

    # need at least 2 individuals/events for rarefaction to be meaningful
    if (n_tot < 2L) next

    # start subsampling from 2 to avoid vegan warnings at size 1
    start <- max(2L, as.integer(step))
    subs  <- seq(start, n_tot, by = step)
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

#' PERMANOVA wrapper (legacy)
#'
#' Maintained for backward compatibility. Delegates to
#' permanova_and_dispersion() for Bray–Curtis vs `type` and returns the
#' tidy table structure used in the practicum.
#'
#' @param dist_obj  stats::dist (typically Bray–Curtis)
#' @param meta_df   data frame with columns site and type
#' @return          tibble with columns distance, term, df, pseudo_F, R2, p_perm, disp_F, disp_p, note
#' @export
permanova_region <- function(dist_obj, meta_df) {
  permanova_and_dispersion(dist_obj, meta_df, group_col = "type", n_perm = 999L, seed = 1L)
}

#' Event counts for focal species across sites
#'
#' @param df           joined data frame with columns site, common, event
#' @param species_vec  character vector of focal species (common values)
#' @return             wide table: site × species, values = event counts
#' @export
species_trend <- function(df, species_vec) {
  df |>
    filter(common %in% species_vec) |>
    count(site, common, name = "events") |>
    pivot_wider(names_from = common, values_from = events, values_fill = 0)
}
