# R/helpers_plots.R

import::here(ggplot, aes, geom_tile, scale_fill_gradient, geom_point,
             geom_text, facet_wrap, theme_minimal, labs, .from = "ggplot2")
import::here(melt, .from = "reshape2")
import::here(select, left_join, filter, .from = "dplyr")
import::here(metaMDS, scores, .from = "vegan")

# Heatmap of Bray–Curtis dissimilarity
plot_turnover_heatmap <- function(dist_obj, meta_df) {
  # distance matrix → long form
  mlt <- melt(as.matrix(dist_obj),
              varnames = c("site1", "site2"),
              value.name = "bray")

  # add region for both sites
  meta_small <- select(meta_df, site, type)

  mlt <- left_join(mlt, meta_small,
                   by = c("site1" = "site")) |>
         left_join(meta_small,
                   by = c("site2" = "site"),
                   suffix = c("_1", "_2"))

  # keep only within-region comparisons
  mlt <- filter(mlt, type_1 == type_2)

  ggplot(mlt, aes(x = site1, y = site2, fill = bray)) +
    geom_tile() +
    facet_wrap(~ type_1, scales = "free", ncol = 1) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Bray–Curtis") +
    theme_minimal() +
    labs(x = "Site", y = "Site")
}

# NMDS ordination
plot_nmds <- function(dist_obj, meta_df, k = 2) {
  ord  <- metaMDS(dist_obj, k = k, trymax = 50, autotransform = FALSE)
  sco  <- as.data.frame(scores(ord, display = "sites"))
  sco$site <- rownames(sco)

  dat  <- left_join(sco, meta_df, by = "site")

  ggplot(dat, aes(x = NMDS1, y = NMDS2, colour = type)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = site), vjust = -0.4, size = 3) +
    theme_minimal() +
    labs(title = "NMDS ordination (Bray–Curtis)",
         colour = "Type")
}
