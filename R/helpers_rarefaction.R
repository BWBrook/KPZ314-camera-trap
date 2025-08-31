# R/helpers_rarefaction.R — rarefaction plots with coverage

import::here(rarefy, estimateR, .from = "vegan")
import::here(ggplot, aes, geom_line, geom_point, theme_minimal, labs,
             scale_colour_gradient, scale_colour_manual, guides, guide_colourbar,
             .from = "ggplot2")
import::here(bind_rows, .from = "dplyr")

#' Site-level rarefaction curves with coverage endpoint
#'
#' @param comm_mat matrix rows=sites, cols=species, integer event counts
#' @param step positive integer step for subsampling
#' @return ggplot
#' @export
rarefaction_site <- function(comm_mat, step = 2L) {
  comm_mat <- as.matrix(comm_mat)
  storage.mode(comm_mat) <- "integer"
  curves <- vector("list", nrow(comm_mat))
  pts    <- vector("list", nrow(comm_mat))

  for (i in seq_len(nrow(comm_mat))) {
    counts <- comm_mat[i, ]
    n_tot  <- sum(counts)
    if (n_tot < 2L) next

    subs  <- seq(2L, n_tot, by = step)
    s_hat <- suppressWarnings(rarefy(counts, sample = subs))
    s_obs <- sum(counts > 0)
    ch1   <- as.numeric(estimateR(counts)["S.chao1"]) 
    covg  <- if (is.finite(ch1) && ch1 > 0) pmin(pmax(s_obs / ch1, 0), 1) else NA_real_

    curves[[i]] <- data.frame(
      site  = rownames(comm_mat)[i],
      reads = subs,
      S     = as.numeric(s_hat)
    )
    pts[[i]] <- data.frame(
      site  = rownames(comm_mat)[i],
      reads = n_tot,
      S     = s_obs,
      coverage = covg
    )
  }

  df  <- bind_rows(curves)
  end <- bind_rows(pts)

  ggplot(df, aes(reads, S, group = site)) +
    geom_line(alpha = 0.6) +
    geom_point(data = end, aes(colour = coverage), size = 2) +
    scale_colour_gradient(name = "Coverage",
                          low = "orange", high = "darkgreen", na.value = "grey50") +
    theme_minimal() +
    labs(title = "Site rarefaction with coverage endpoints",
         x = "Number of events (subsample)", y = "Expected richness")
}

#' Group-level pooled rarefaction with coverage annotations
#'
#' @param comm_mat site × species matrix (rows named by site)
#' @param meta_df  data frame mapping site -> group (column group_col)
#' @param group_col grouping column in meta_df (default "type")
#' @param step subsampling step
#' @return ggplot
#' @export
rarefaction_by_group <- function(comm_mat, meta_df, group_col = "type", step = 5L) {
  stopifnot(group_col %in% names(meta_df))
  comm_mat <- as.matrix(comm_mat)
  storage.mode(comm_mat) <- "integer"
  groups <- unique(meta_df[[group_col]])

  curves <- list(); ann <- list()
  for (g in groups) {
    sites_g <- meta_df$site[meta_df[[group_col]] == g]
    idx <- intersect(rownames(comm_mat), as.character(sites_g))
    if (length(idx) == 0L) next
    v <- colSums(comm_mat[idx, , drop = FALSE])
    n_tot <- sum(v)
    if (n_tot < 2L) next
    subs <- seq(2L, n_tot, by = step)
    s_hat <- suppressWarnings(rarefy(v, sample = subs))
    s_obs <- sum(v > 0)
    ch1   <- as.numeric(estimateR(v)["S.chao1"]) 
    covg  <- if (is.finite(ch1) && ch1 > 0) pmin(pmax(s_obs / ch1, 0), 1) else NA_real_
    curves[[as.character(g)]] <- data.frame(group = as.character(g), reads = subs, S = as.numeric(s_hat))
    ann[[as.character(g)]]    <- data.frame(group = as.character(g), reads = n_tot, S = s_obs, coverage = covg)
  }
  df <- bind_rows(curves)
  ad <- bind_rows(ann)

  ggplot(df, aes(reads, S, colour = group)) +
    geom_line(linewidth = 1.1) +
    geom_point(data = ad, aes(fill = coverage), colour = "black", shape = 21, size = 2) +
    scale_colour_manual(values = c(dry = "sienna3", wet = "steelblue4"), name = group_col) +
    scale_fill_gradient(name = "Coverage",
                        low = "orange", high = "darkgreen", na.value = "grey50") +
    theme_minimal() +
    labs(title = "Group rarefaction with coverage endpoints",
         x = "Number of events (subsample)", y = "Expected richness")
}
