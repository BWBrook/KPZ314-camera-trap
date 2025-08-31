# R/helpers_metrics.R

import::here(group_by, summarise, n, n_distinct, ungroup, tibble, mutate, select, left_join, distinct, .from = "dplyr")
import::here(pivot_wider, replace_na, .from = "tidyr")
import::here(diversity, estimateR, .from = "vegan")
import::here(abort, .from = "rlang")

#' Alpha‑diversity metrics per site
#'
#' @param df Data frame containing at minimum:
#'   * site  (factor or character)
#'   * common   (factor or character)
#'   * event        (integer – independent events)
#'
#' @return tibble with one row per site and columns:
#'   richness, shannon, simpson, pielou, chao1, hill_q0, hill_q1, hill_q2
#'
#' @export
calc_alpha <- function(df) {
  ## collapse to a site × species abundance table (event counts)
  mat <- df |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site, common) |>
    summarise(events = n(), .groups = "drop") |>
    pivot_wider(names_from = common,
                values_from = events,
                values_fill = 0L)

  site_id <- mat$site
  mat_num <- as.matrix(mat[ , -1, drop = FALSE ])

  richness <- rowSums(mat_num > 0)
  shannon  <- diversity(mat_num, index = "shannon")
  simpson  <- diversity(mat_num, index = "simpson")
  pielou   <- shannon / log(pmax(richness, 1))
  chao1    <- unname(apply(mat_num, 1, function(v) estimateR(v)["S.chao1"]))
  hill_q0  <- richness
  hill_q1  <- exp(shannon)
  hill_q2  <- diversity(mat_num, index = "invsimpson")

  tibble(
    site = site_id,
    richness    = richness,
    shannon     = shannon,
    simpson     = simpson,
    pielou      = pielou,
    chao1       = chao1,
    hill_q0     = hill_q0,
    hill_q1     = hill_q1,
    hill_q2     = hill_q2
  )
}

#' Bootstrap CIs for alpha diversity (per site)
#'
#' Resamples events within each site to obtain uncertainty on Hill numbers
#' and Pielou's evenness.
#'
#' @param df_events Data frame with site, common, event
#' @param n_boot Number of bootstrap replicates (default 500)
#' @param seed RNG seed for reproducibility
#' @return tibble(site, metric, est, lwr, upr, se, n_boot)
#' @export
bootstrap_alpha <- function(df_events, n_boot = 500L, seed = 1L) {
  if (!all(c("site","common","event") %in% names(df_events))) {
    abort("bootstrap_alpha(): df_events must contain site, common, event.")
  }

  # per-site species event counts
  wide <- df_events |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site, common) |>
    summarise(events = n(), .groups = "drop")

  sites <- unique(wide$site)
  if (length(sites) == 0L) abort("bootstrap_alpha(): no sites found.")

  # point estimates
  est_tbl <- calc_alpha(df_events) |>
    select(site, hill_q0, hill_q1, hill_q2, pielou)

  set.seed(as.integer(seed))
  res_list <- lapply(split(wide, wide$site), function(df_site) {
    counts <- df_site$events
    sp     <- as.character(df_site$common)
    names(counts) <- sp
    total <- sum(counts)
    if (total <= 0) abort("bootstrap_alpha(): site has zero events.")
    labels <- rep(sp, times = counts)

    boot_vals <- replicate(n_boot, {
      samp  <- sample(labels, size = length(labels), replace = TRUE)
      tab   <- table(samp)
      v     <- as.numeric(tab)
      names(v) <- names(tab)
      # richness
      richness <- sum(v > 0)
      # shannon, invsimpson from vegan::diversity
      sh <- diversity(v, index = "shannon")
      q1 <- exp(sh)
      q2 <- diversity(v, index = "invsimpson")
      pielou <- if (richness > 0) sh / log(richness) else NA_real_
      c(hill_q0 = richness, hill_q1 = q1, hill_q2 = q2, pielou = pielou)
    })

    metrics <- rownames(boot_vals)
    lwr <- apply(boot_vals, 1, function(x) quantile(x, 0.025, na.rm = TRUE))
    est <- apply(boot_vals, 1, function(x) mean(x, na.rm = TRUE))
    upr <- apply(boot_vals, 1, function(x) quantile(x, 0.975, na.rm = TRUE))
    se  <- apply(boot_vals, 1, sd, na.rm = TRUE)

    tibble(
      site = df_site$site[1],
      metric = metrics,
      lwr = as.numeric(lwr), est = as.numeric(est), upr = as.numeric(upr), se = as.numeric(se),
      n_boot = as.integer(n_boot)
    )
  })

  boot_tbl <- do.call(rbind, res_list)

  # replace est with point estimates from calc_alpha
  est_long <- est_tbl |>
    tidyr::pivot_longer(-site, names_to = "metric", values_to = "est_point")

  left_join(boot_tbl, est_long, by = c("site","metric")) |>
    mutate(est = est_point) |>
    select(-est_point)
}

#' Gamma Hill numbers with site-bootstrap CIs
#'
#' @param df_events Data frame with site, common, event and (optionally) 'type'
#' @param group "all" for overall, or a column name present in df (e.g., "type")
#' @param qs vector of q orders (supports 0,1,2)
#' @param n_boot bootstrap replicates
#' @param seed RNG seed
#' @return tibble(group, level, q, est, lwr, upr, se, s_obs, chao1, coverage)
#' @export
gamma_hill <- function(df_events, group = c("all","type"), qs = c(0L,1L,2L),
                       n_boot = 1000L, seed = 1L) {
  group <- match.arg(group)
  if (!all(c("site","common","event") %in% names(df_events))) {
    abort("gamma_hill(): df_events must contain site, common, event.")
  }

  # collapse to site × species counts
  counts_df <- df_events |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site, common) |>
    summarise(events = n(), .groups = "drop")

  set.seed(as.integer(seed))

  do_group <- function(level_label, site_ids) {
    # observed pooled counts across sites in this level
    pool <- counts_df[counts_df$site %in% site_ids, ] |>
      group_by(common) |>
      summarise(events = sum(events), .groups = "drop")
    v <- pool$events
    names(v) <- as.character(pool$common)

    s_obs  <- sum(v > 0)
    ch1    <- as.numeric(estimateR(v)["S.chao1"]) 
    covg   <- if (is.finite(ch1) && ch1 > 0) pmin(pmax(s_obs / ch1, 0), 1) else NA_real_

    # function to compute hill numbers from counts
    hill_from_counts <- function(v, q) {
      if (q == 0) return(sum(v > 0))
      if (q == 1) return(exp(diversity(v, index = "shannon")))
      if (q == 2) return(diversity(v, index = "invsimpson"))
      NA_real_
    }

    # site-level bootstrap: resample sites with replacement and pool counts
    sites_level <- unique(site_ids)
    if (length(sites_level) == 0L) abort("gamma_hill(): group has zero sites.")

    boot_vals <- lapply(qs, function(q) {
      replicate(n_boot, {
        samp_sites <- sample(sites_level, size = length(sites_level), replace = TRUE)
        pool_b <- counts_df[counts_df$site %in% samp_sites, ] |>
          group_by(common) |>
          summarise(events = sum(events), .groups = "drop")
        vv <- pool_b$events
        names(vv) <- as.character(pool_b$common)
        hill_from_counts(vv, q)
      })
    })

    ests <- vapply(qs, function(q) hill_from_counts(v, q), numeric(1))
    tibble(
      q   = qs,
      est = ests,
      lwr = vapply(boot_vals, function(x) quantile(x, 0.025, na.rm = TRUE), numeric(1)),
      upr = vapply(boot_vals, function(x) quantile(x, 0.975, na.rm = TRUE), numeric(1)),
      se  = vapply(boot_vals, function(x) sd(x, na.rm = TRUE), numeric(1)),
      s_obs = s_obs,
      chao1 = ch1,
      coverage = covg
    ) |>
      mutate(level = level_label)
  }

  if (group == "all") {
    site_ids <- unique(counts_df$site)
    out <- do_group("overall", site_ids)
    out$group <- "all"
    return(out[, c("group","level","q","est","lwr","upr","se","s_obs","chao1","coverage")])
  }

  if (!("type" %in% names(df_events))) {
    abort("gamma_hill(): 'type' column required for group = 'type'.")
  }
  st <- distinct(df_events, site, type)
  sites_by_level <- split(st$site, st$type)
  res <- lapply(names(sites_by_level), function(lvl) do_group(lvl, sites_by_level[[lvl]]))
  out <- do.call(rbind, res)
  out$group <- "type"
  out[, c("group","level","q","est","lwr","upr","se","s_obs","chao1","coverage")]
}

#' Community matrix (sites × species, abundance = event count)
#'
#' @param df Data frame with site, common, event
#' @return matrix; row names = site, col names = common
#' @export
build_comm_matrix <- function(df) {
  wide <- df |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site, common) |>
    summarise(events = n(), .groups = "drop") |>
    pivot_wider(
      names_from  = common,
      values_from = events,
      values_fill = 0L
    )

  mat <- as.matrix(wide[ , -1, drop = FALSE ])
  rownames(mat) <- wide$site
  mat
}
