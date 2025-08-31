# R/helpers_plots.R

import::here(ggplot, aes, geom_tile, scale_fill_gradient, geom_point,
             geom_text, facet_wrap, theme_minimal, labs, coord_sf,
             scale_colour_manual, coord_quickmap, geom_sf, .from = "ggplot2")
import::here(melt, .from = "reshape2")
import::here(select, left_join, filter, .from = "dplyr")
import::here(metaMDS, scores, .from = "vegan")
import::here(st_as_sf, st_transform, st_bbox, st_as_sfc, st_buffer, .from = "sf")

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

# Site map with optional OSM basemap (if ggspatial is installed)
plot_site_map <- function(meta_df, buffer_m = 500L) {
  stopifnot(all(c("lat", "lon", "type") %in% names(meta_df)))

  has_tiles <- "ggspatial" %in% rownames(utils::installed.packages())

  if (has_tiles) {
    # Work in Web Mercator for tiles, compute buffered extent in meters
    pts_ll <- st_as_sf(meta_df, coords = c("lon", "lat"), crs = 4326)
    pts_m  <- st_transform(pts_ll, 3857)
    bb_m   <- st_bbox(st_buffer(st_as_sfc(st_bbox(pts_m)), dist = as.numeric(buffer_m)))

    ggplot() +
      ggspatial::annotation_map_tile(zoomin = 0, progress = "none") +
      geom_sf(data = pts_m, aes(colour = type), size = 3) +
      coord_sf(crs = 3857,
               xlim = c(bb_m[["xmin"]], bb_m[["xmax"]]),
               ylim = c(bb_m[["ymin"]], bb_m[["ymax"]]),
               expand = FALSE) +
      theme_minimal() +
      labs(title = "Camera sites") +
      scale_colour_manual(values = c(dry = "sienna3", wet = "steelblue4"),
                          name = "Type")
  } else {
    # Fallback: lon/lat scatter with equal-aspect mapping and buffered bounds
    mean_lat <- mean(meta_df$lat, na.rm = TRUE)
    lat_buf_deg <- buffer_m / 111320
    lon_buf_deg <- buffer_m / (111320 * cos(mean_lat * pi / 180))
    xlim <- range(meta_df$lon, na.rm = TRUE) + c(-lon_buf_deg, lon_buf_deg)
    ylim <- range(meta_df$lat, na.rm = TRUE) + c(-lat_buf_deg, lat_buf_deg)

    ggplot(meta_df, aes(x = lon, y = lat, colour = type)) +
      geom_point(size = 3) +
      coord_quickmap(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_minimal() +
      labs(title = "Camera sites", x = "Longitude", y = "Latitude") +
      scale_colour_manual(values = c(dry = "sienna3", wet = "steelblue4"),
                          name = "Type")
  }
}
