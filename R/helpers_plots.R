# R/helpers_plots.R

import::here(ggplot, aes, geom_tile, scale_fill_gradient, geom_point,
             geom_text, facet_wrap, theme_minimal, labs, coord_sf,
             scale_colour_manual, coord_quickmap, geom_sf, geom_polygon,
             guides, guide_none, geom_pointrange, scale_x_log10, geom_vline,
             .from = "ggplot2")
import::here(melt, .from = "reshape2")
import::here(select, left_join, filter, tibble, .from = "dplyr")
import::here(metaMDS, scores, .from = "vegan")
import::here(ordiellipse, .from = "vegan")
import::here(st_as_sf, st_transform, st_bbox, st_as_sfc, st_buffer, .from = "sf")
import::here(local_seed, .from = "withr")
import::here(cli_warn, .from = "cli")
import::here(tibble, .from = "dplyr")
import::here(melt, .from = "reshape2")

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

# General heatmap by within-group facets (works for any distance)
plot_heat_by_group <- function(dist_obj, meta_df, group_col = "type") {
  mat <- as.matrix(dist_obj)
  # Build a global ordering: within each group, order sites by hclust on the submatrix
  meta_small0 <- meta_df[, c("site", group_col)]
  names(meta_small0) <- c("site", "group")
  ord_levels <- character(0)
  for (g in unique(meta_small0$group)) {
    sites_g <- meta_small0$site[meta_small0$group == g]
    sites_g <- intersect(sites_g, rownames(mat))
    if (length(sites_g) >= 2L) {
      sub <- as.dist(mat[sites_g, sites_g, drop = FALSE])
      hc <- stats::hclust(sub, method = "average")
      ord_levels <- c(ord_levels, sites_g[hc$order])
    } else {
      ord_levels <- c(ord_levels, sites_g)
    }
  }

  mlt <- melt(mat,
              varnames = c("site1", "site2"),
              value.name = "distance")

  # join group for both sites (base indexing to avoid NSE)
  if (!("site" %in% names(meta_df))) rlang::abort("plot_heat_by_group(): meta_df must contain 'site'.")
  if (!(group_col %in% names(meta_df))) rlang::abort(sprintf("plot_heat_by_group(): meta_df must contain '%s'.", group_col))
  meta_small <- meta_df[, c("site", group_col)]
  names(meta_small) <- c("site", "group")

  mlt <- left_join(mlt, meta_small, by = c("site1" = "site")) |>
    left_join(meta_small, by = c("site2" = "site"), suffix = c("_1", "_2"))

  # warn if any group has <2 sites
  counts <- table(meta_small$group)
  if (any(counts < 2L)) cli_warn("One or more groups have <2 sites; panels may be sparse.")

  mlt <- filter(mlt, group_1 == group_2)
  # apply global ordering for readability
  mlt$site1 <- factor(mlt$site1, levels = ord_levels)
  mlt$site2 <- factor(mlt$site2, levels = ord_levels)
  method <- attr(dist_obj, "method"); if (is.null(method)) method <- "distance"
  method_chr <- as.character(method)
  cap <- if (method_chr %in% c("bray","jaccard")) turnover_caption(method_chr) else turnover_caption()

  ggplot(mlt, aes(x = site1, y = site2, fill = distance)) +
    geom_tile() +
    facet_wrap(~ group_1, scales = "free", ncol = 1) +
    scale_fill_gradient(low = "white", high = "red", name = method) +
    theme_minimal() +
    labs(x = "Site", y = "Site", caption = cap)
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

# NMDS with convex hulls and stats (seeded)
nmds_with_hulls <- function(dist_obj, meta_df, group_col = "type", k = 2L, seed = 1L) {
  if (!("site" %in% names(meta_df))) rlang::abort("nmds_with_hulls(): meta_df must contain 'site'.")
  if (!(group_col %in% names(meta_df))) rlang::abort(sprintf("nmds_with_hulls(): meta_df must contain '%s'.", group_col))

  local_seed(as.integer(seed))
  ord <- try(metaMDS(dist_obj, k = as.integer(k), trymax = 50, autotransform = FALSE), silent = TRUE)
  if (inherits(ord, "try-error")) {
    rlang::abort("nmds_with_hulls(): NMDS failed; consider k = 3.")
  }
  sco <- as.data.frame(scores(ord, display = "sites"))
  sco$site <- rownames(sco)

  meta_small <- meta_df[, c("site", group_col)]
  names(meta_small) <- c("site", "group")
  dat <- left_join(sco, meta_small, by = "site")

  # build convex hulls per group when >= 3 points
  hulls <- lapply(split(dat, dat$group), function(df) {
    if (nrow(df) < 3L) return(NULL)
    idx <- grDevices::chull(df$NMDS1, df$NMDS2)
    df[idx, c("NMDS1","NMDS2","group"), drop = FALSE]
  })
  hulls <- do.call(rbind, hulls)

  method_chr <- as.character(attr(dist_obj, "method"))
  cap <- if (method_chr %in% c("bray","jaccard")) turnover_caption(method_chr) else turnover_caption()

  p <- ggplot(dat, aes(x = NMDS1, y = NMDS2, colour = group)) +
    { if (!is.null(hulls)) geom_polygon(
        data = hulls,
        aes(x = NMDS1, y = NMDS2, group = group, fill = group),
        colour = NA, alpha = 0.2, inherit.aes = FALSE
      ) else NULL } +
    geom_point(size = 3, alpha = 0.85) +
    geom_text(aes(label = site), vjust = -0.4, size = 3, show.legend = FALSE) +
    theme_minimal() +
    labs(title = paste0("NMDS (", method_chr, ") with convex hulls"), colour = group_col, caption = cap)

  st <- tibble(
    method = as.character(attr(dist_obj, "method")),
    k = as.integer(k),
    stress = as.numeric(ord$stress),
    warn = ifelse(as.numeric(ord$stress) > 0.2, "stress>0.2: interpret cautiously", ""),
    seed = as.integer(seed)
  )

  list(fig = p, stats = st)
}

# NMDS with ellipses (mean ± 1 SD) as an alternative to hulls
nmds_with_ellipses <- function(dist_obj, meta_df, group_col = "type", k = 2L, seed = 1L) {
  if (!("site" %in% names(meta_df))) rlang::abort("nmds_with_ellipses(): meta_df must contain 'site'.")
  if (!(group_col %in% names(meta_df))) rlang::abort(sprintf("nmds_with_ellipses(): meta_df must contain '%s'.", group_col))
  local_seed(as.integer(seed))
  ord <- try(metaMDS(dist_obj, k = as.integer(k), trymax = 50, autotransform = FALSE), silent = TRUE)
  if (inherits(ord, "try-error")) rlang::abort("nmds_with_ellipses(): NMDS failed; consider k = 3.")
  sco <- as.data.frame(scores(ord, display = "sites"))
  sco$site <- rownames(sco)
  meta_small <- meta_df[, c("site", group_col)]
  names(meta_small) <- c("site", "group")
  dat <- left_join(sco, meta_small, by = "site")
  method_chr <- as.character(attr(dist_obj, "method"))
  cap <- if (method_chr %in% c("bray","jaccard")) turnover_caption(method_chr) else turnover_caption()
  p <- ggplot(dat, aes(NMDS1, NMDS2, colour = group)) +
    geom_point(size = 3, alpha = 0.85) +
    theme_minimal() +
    labs(title = paste0("NMDS (", method_chr, ") with ordiellipse"), colour = group_col, caption = cap)
  # simple ellipse via ggplot (normal approximation)
  p <- p + ggplot2::stat_ellipse(aes(group = group, colour = group), type = "norm")
  list(fig = p, stats = tibble(method = method_chr, k = as.integer(k), stress = as.numeric(ord$stress), warn = ifelse(as.numeric(ord$stress) > 0.2, "stress>0.2", ""), seed = as.integer(seed)))
}
# Forest plot for GLM IRRs (non-intercept terms)
plot_glm_coefs <- function(irrs_tbl) {
  df <- irrs_tbl
  if (!all(c("term","irr","lwr","upr","species") %in% names(df)))
    rlang::abort("plot_glm_coefs(): missing required columns in irrs_tbl.")

  df <- df[df$term != "(Intercept)", , drop = FALSE]
  # readable term labels
  lab_map <- c(
    typewet = "wet vs dry",
    leaf_litter_z = "leaf litter (+1 SD)",
    cwd_z = "CWD (+1 SD)",
    can_foliage_z = "canopy foliage (+1 SD)"
  )
  df$label <- lab_map[df$term]
  df$label[is.na(df$label)] <- df$term

  ggplot(df, aes(x = irr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = lwr, xmax = upr)) +
    scale_x_log10() +
    theme_minimal() +
    labs(
      title = as.character(df$species[1]),
      x = "Incidence Rate Ratio (IRR)", y = NULL,
      subtitle = "Numeric covariates are per +1 SD"
    )
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

# Detection history heatmap (wrapper retained for consistency)
# Prefer using plot_detection_heatmap() defined in helpers_detectability.R.

# psi curves plot helper (wrapper); the core implementation lives in helpers_detectability.R
