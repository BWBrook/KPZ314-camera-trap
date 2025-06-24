# R/helpers_plots.R
# Heatmap of Bray–Curtis dissimilarity
plot_turnover_heatmap <- function(dist_obj, meta_df) {
  import::from("ggplot2", ggplot, aes, geom_tile, scale_fill_gradient,
               facet_wrap, theme_minimal, labs)
  import::from("reshape2", melt)
  import::from("dplyr", select, left_join, filter)

  # distance matrix → long form
  mlt <- melt(as.matrix(dist_obj),
              varnames = c("site1", "site2"),
              value.name = "bray")

  # add region for both sites
  meta_small <- select(meta_df, camera_site, region)

  mlt <- left_join(mlt, meta_small,
                   by = c("site1" = "camera_site")) |>
         left_join(meta_small,
                   by = c("site2" = "camera_site"),
                   suffix = c("_1", "_2"))

  # keep only within-region comparisons
  mlt <- filter(mlt, region_1 == region_2)

  ggplot(mlt, aes(x = site1, y = site2, fill = bray)) +
    geom_tile() +
    facet_wrap(~ region_1, scales = "free") +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Bray–Curtis") +
    theme_minimal() +
    labs(x = "Site", y = "Site")
}

# NMDS ordination
plot_nmds <- function(dist_obj, meta_df, k = 2) {
  import::from("vegan", metaMDS, scores)
  import::from("ggplot2", ggplot, aes, geom_point, geom_text, theme_minimal,
               labs)
  import::from("dplyr", left_join)

  ord  <- metaMDS(dist_obj, k = k, trymax = 50, autotransform = FALSE)
  sco  <- as.data.frame(scores(ord, display = "sites"))
  sco$camera_site <- rownames(sco)

  dat  <- left_join(sco, meta_df, by = "camera_site")

  ggplot(dat, aes(x = NMDS1, y = NMDS2, colour = region)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(aes(label = camera_site), vjust = -0.4, size = 3) +
    theme_minimal() +
    labs(title = "NMDS ordination (Bray–Curtis)",
         colour = "Region")
}
