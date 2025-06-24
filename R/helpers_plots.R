# R/helpers_plots.R
# Heatmap of Bray–Curtis dissimilarity
plot_turnover_heatmap <- function(dist_obj, meta_df) {
  import::from("ggplot2", ggplot, aes, geom_tile, scale_fill_gradient,
               facet_wrap, theme_minimal, labs)
  import::from("reshape2", melt)

  mat   <- as.matrix(dist_obj)
  mlt   <- melt(mat, varnames = c("site1", "site2"))
  mlt   <- merge(mlt, meta_df[ , c("camera_site", "region") ],
                 by.x = "site1", by.y = "camera_site", all.x = TRUE)

  ggplot(mlt, aes(x = site1, y = site2, fill = value)) +
    geom_tile() +
    facet_wrap(~ region, scales = "free") +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    labs(fill = "Bray–Curtis",
         x = "Site",
         y = "Site")
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
