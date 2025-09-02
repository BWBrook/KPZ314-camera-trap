# R/helpers_detectability.R — detection histories, naïve occupancy, and sensitivity

import::here(
  mutate, select, left_join, group_by, summarise, n, arrange, ungroup,
  distinct, filter, tibble, bind_rows, .from = "dplyr"
)
import::here(pivot_longer, pivot_wider, .from = "tidyr")
import::here(tibble, .from = "tibble")
import::here(ggplot, aes, geom_tile, geom_line, geom_point, geom_text,
             geom_errorbar, facet_wrap, theme_minimal, labs, scale_y_continuous,
             scale_x_continuous, annotate, guides, guide_none, .from = "ggplot2")
import::here(melt, .from = "reshape2")
import::here(abort, .from = "rlang")
import::from("lubridate", with_tz)
# Explicit imports to make dependency visible to renv; we still call via unmarked:: to avoid generic clashes
import::from("unmarked", occu, unmarkedFrameOccu, backTransform)

# 1) daily detection history for a species
detection_history_for_species <- function(events_df, site_df, species, max_days) {
  req_ev <- c("site","common","datetime")
  if (!all(req_ev %in% names(events_df))) abort("detection_history_for_species(): events_df must have site, common, datetime.")
  req_meta <- c("site","type","first_image","op_days")
  if (!all(req_meta %in% names(site_df))) abort("detection_history_for_species(): site_df must have site, type, first_image, op_days.")

  max_days <- as.integer(max_days)
  if (!is.finite(max_days) || max_days < 1L) abort("detection_history_for_species(): max_days must be >= 1.")

  # Dates sanity
  bad_sites <- site_df$site[is.na(site_df$first_image)]
  if (length(bad_sites)) abort(paste0("detection_history_for_species(): missing first_image for sites: ", paste(bad_sites, collapse = ", ")))

  # Build base matrix
  sites <- as.character(site_df$site)
  n_sites <- length(sites)
  mat <- matrix(0L, nrow = n_sites, ncol = max_days,
                dimnames = list(sites, paste0("d", seq_len(max_days))))

  # Filter species
  ev <- events_df[events_df$common == species, c("site","datetime"), drop = FALSE]
  if (nrow(ev) > 0) {
    # Join first_image to compute day index
    meta_small <- site_df[, c("site","first_image"), drop = FALSE]
    # robust date conversion for first_image in Hobart local time
    first_posix <- try(as.POSIXct(meta_small$first_image, format = "%d/%m/%Y %H:%M:%S", tz = "UTC"), silent = TRUE)
    if (inherits(first_posix, "try-error") || all(is.na(first_posix))) {
      # fall back: treat as Date already
      first_posix <- as.POSIXct(meta_small$first_image, tz = "UTC")
    }
    first_hobart_date <- as.Date(with_tz(first_posix, tzone = "Australia/Hobart"))
    meta_small$first_date <- first_hobart_date
    ev <- left_join(ev, meta_small[, c("site","first_date")], by = "site")
    # day index: 1..max_days using Hobart-local date boundaries
    ev_day_date <- as.Date(with_tz(ev$datetime, tzone = "Australia/Hobart"))
    ev$day <- as.integer(ev_day_date - ev$first_date) + 1L
    ev <- ev[is.finite(ev$day) & ev$day >= 1L & ev$day <= max_days, , drop = FALSE]
    if (nrow(ev) > 0) {
      # any detection on day -> mark 1
      by_site_day <- unique(ev[, c("site","day")])
      for (i in seq_len(nrow(by_site_day))) {
        s <- as.character(by_site_day$site[i])
        d <- as.integer(by_site_day$day[i])
        if (!is.na(match(s, sites))) mat[s, d] <- 1L
      }
    }
  }

  # carry species as attribute for downstream functions
  attr(mat, "species") <- species

  list(
    species = as.character(species),
    days    = seq_len(max_days),
    mat     = mat,
    meta    = tibble(site = sites, type = as.character(site_df$type))
  )
}

# 2) cumulative detection curves by habitat (proportion of sites detected by day t)
cumulative_detection_curves <- function(hist_mat, meta_df) {
  if (!is.matrix(hist_mat)) abort("cumulative_detection_curves(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("cumulative_detection_curves(): meta_df must have site and type.")
  species <- attr(hist_mat, "species") %||% "(unknown)"

  # align rows of hist to meta
  meta_use <- distinct(meta_df[, c("site","type")])
  rn <- rownames(hist_mat)
  if (is.null(rn)) rn <- meta_use$site
  meta_use <- meta_use[match(rn, meta_use$site), , drop = FALSE]
  if (anyNA(meta_use$site)) abort("cumulative_detection_curves(): rownames(hist_mat) must match meta_df$site.")

  max_days <- ncol(hist_mat)
  out <- lapply(seq_len(max_days), function(t) {
    det <- rowSums(hist_mat[, seq_len(t), drop = FALSE]) > 0
    tibble(
      species = species,
      type    = as.character(meta_use$type),
      day     = as.integer(t),
      detected = as.integer(det)
    ) |>
      group_by(species, type, day) |>
      summarise(psi_naive = mean(detected), .groups = "drop")
  })
  do.call(rbind, out)
}

# 3) alias that enforces monotonicity across days
naive_psi_curves <- function(hist_mat, meta_df) {
  tbl <- cumulative_detection_curves(hist_mat, meta_df)
  # monotone non-decreasing by species x type
  chk <- tbl |>
    group_by(species, type) |>
    arrange(day, .by_group = TRUE) |>
    summarise(non_dec = all(diff(psi_naive) >= -1e-12), .groups = "drop")
  if (!all(chk$non_dec)) abort("naive_psi_curves(): psi(t) must be non-decreasing; check detection history logic.")
  tbl
}

# 4) heatmap of detection history (0/1) faceted by habitat
plot_detection_heatmap <- function(hist_mat, meta_df, species) {
  if (!is.matrix(hist_mat)) abort("plot_detection_heatmap(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("plot_detection_heatmap(): meta_df must have site and type.")

  # long form
  df <- melt(hist_mat, varnames = c("site","day"), value.name = "det")
  # convert day to integer index 1..T
  df$day <- as.integer(sub("^d", "", as.character(df$day)))
  meta_small <- distinct(meta_df[, c("site","type")])
  df <- left_join(df, meta_small, by = "site")

  # ordering by first detection day (Inf if never)
  pos <- df[df$det > 0, , drop = FALSE]
  if (nrow(pos) > 0) {
    fd <- stats::aggregate(day ~ site, data = pos, FUN = min)
    fd <- left_join(tibble(site = unique(df$site)), fd, by = "site")
    fd$day[is.na(fd$day)] <- Inf
  } else {
    fd <- tibble(site = unique(df$site), day = Inf)
  }
  ord <- df |>
    distinct(site, type) |>
    left_join(fd, by = "site") |>
    arrange(type, day, site)
  df$site <- factor(df$site, levels = ord$site)

  p <- ggplot(df, aes(x = day, y = site, fill = factor(det))) +
    geom_tile() +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = pretty(seq_len(ncol(hist_mat)))) +
    ggplot2::scale_y_discrete(NULL, breaks = NULL) +
    theme_minimal() +
    labs(
      title = paste0("Daily detections — ", species),
      x = "Day since first_image",
      caption = "Days without detections treated as non-detections; cameras assumed operating continuously for op_days."
    ) +
    guides(fill = "none")
  attr(p, "species") <- species
  p
}

# 5) psi-hat(t) curves with 0.8 threshold annotations
plot_psi_curves <- function(curves_tbl) {
  req <- c("species","type","day","psi_naive")
  if (!all(req %in% names(curves_tbl))) abort("plot_psi_curves(): curves_tbl missing required columns.")
  sp <- as.character(unique(curves_tbl$species))
  p <- ggplot(curves_tbl, aes(x = day, y = psi_naive, colour = type)) +
    geom_line(linewidth = 1) +
    theme_minimal() +
    scale_y_continuous(limits = c(0,1)) +
    labs(title = paste0("Naïve occupancy ψ̂(t) — ", sp), x = "Days", y = "ψ̂(t)")

  # annotate first day reaching 0.8 per type (if any)
  ann <- curves_tbl |>
    group_by(type) |>
    arrange(day, .by_group = TRUE) |>
    summarise(d80 = if (any(psi_naive >= 0.8)) min(day[psi_naive >= 0.8]) else NA_integer_, .groups = "drop")
  if (nrow(ann)) {
    ann_use <- ann[is.finite(ann$d80) & !is.na(ann$d80), , drop = FALSE]
    if (nrow(ann_use)) {
      p <- p + ggplot2::geom_vline(xintercept = ann_use$d80, linetype = "dotted", alpha = 0.6)
      # label text near top
      lab_df <- data.frame(
        day = ann_use$d80,
        psi_naive = 0.95,
        type = ann_use$type,
        lab = "~80% of occupied sites detected (p<1)"
      )
      p <- p + geom_text(data = lab_df, aes(x = day, y = psi_naive, label = lab, colour = type),
                         inherit.aes = FALSE, size = 3, vjust = -0.3, show.legend = FALSE)
    }
  }
  p <- p + labs(caption = "Dashed line marks first day ψ̂(t) ≥ 0.8; this is not '80% occupancy'.")
  p
}

# 6) remove unknown_animal events
events_without_unknown <- function(events_df) {
  if (!("common" %in% names(events_df))) abort("events_without_unknown(): events_df must contain 'common'.")
  events_df[events_df$common != "unknown_animal", , drop = FALSE]
}

# 7) alpha deltas with and without unknown_animal
alpha_delta_unknown <- function(alpha_all, alpha_no_unk) {
  metrics <- c("richness","shannon","simpson","pielou","chao1","hill_q1","hill_q2")
  if (!all(c("site", metrics) %in% names(alpha_all))) abort("alpha_delta_unknown(): alpha_all missing required columns.")
  if (!all(c("site", metrics) %in% names(alpha_no_unk))) abort("alpha_delta_unknown(): alpha_no_unk missing required columns.")

  long_a <- pivot_longer(alpha_all[, c("site", metrics)], -site, names_to = "metric", values_to = "all")
  long_b <- pivot_longer(alpha_no_unk[, c("site", metrics)], -site, names_to = "metric", values_to = "no_unk")
  out <- left_join(long_a, long_b, by = c("site","metric")) |>
    mutate(delta = no_unk - all) |>
    select(site, metric, delta)

  # richness must not increase when removing unknowns
  bad <- out[out$metric == "richness" & out$delta > 1e-12, , drop = FALSE]
  if (nrow(bad)) abort("alpha_delta_unknown(): richness delta > 0 for some sites; check duplicate handling.")
  out
}

# 8) beta summary deltas comparing all vs no-unknown for Bray and Jaccard
beta_delta_summary <- function(comm_all, comm_no_unk) {
  import::from("stats", "cor")
  # helper to compute upper triangle vector
  ut <- function(m) as.numeric(m[upper.tri(m, diag = FALSE)])

  # Bray
  db_all <- compute_dist(comm_all, method = "bray")
  db_nu  <- compute_dist(comm_no_unk, method = "bray")
  mb_all <- as.matrix(db_all); mb_nu <- as.matrix(db_nu)
  # align by shared site labels
  labs <- intersect(rownames(mb_all), rownames(mb_nu))
  mb_all <- mb_all[labs, labs, drop = FALSE]
  mb_nu  <- mb_nu[labs, labs, drop = FALSE]
  vb_all <- ut(mb_all); vb_nu <- ut(mb_nu)
  r_b <- suppressWarnings(cor(vb_all, vb_nu, use = "complete.obs"))
  row_b <- tibble(distance = "bray",
                  r_upper = as.numeric(r_b),
                  mean_abs_diff = mean(abs(vb_all - vb_nu)),
                  median_abs_diff = stats::median(abs(vb_all - vb_nu)),
                  n_pairs = length(vb_all))

  # Jaccard
  dj_all <- compute_dist(comm_all, method = "jaccard")
  dj_nu  <- compute_dist(comm_no_unk, method = "jaccard")
  mj_all <- as.matrix(dj_all); mj_nu <- as.matrix(dj_nu)
  labs2 <- intersect(rownames(mj_all), rownames(mj_nu))
  mj_all <- mj_all[labs2, labs2, drop = FALSE]
  mj_nu  <- mj_nu[labs2, labs2, drop = FALSE]
  vj_all <- ut(mj_all); vj_nu <- ut(mj_nu)
  r_j <- suppressWarnings(cor(vj_all, vj_nu, use = "complete.obs"))
  row_j <- tibble(distance = "jaccard",
                  r_upper = as.numeric(r_j),
                  mean_abs_diff = mean(abs(vj_all - vj_nu)),
                  median_abs_diff = stats::median(abs(vj_all - vj_nu)),
                  n_pairs = length(vj_all))

  bind_rows(row_b, row_j)
}

# 9) unmarked frame builder
occu_build_umf <- function(hist_mat, meta_df) {
  if (!is.matrix(hist_mat)) abort("occu_build_umf(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("occu_build_umf(): meta_df must have site and type.")
  if (!requireNamespace("unmarked", quietly = TRUE)) {
    abort("Package 'unmarked' is required for the occupancy demo. Add it via dependencies.R then renv::restore().")
  }
  # enforce row alignment
  rn <- rownames(hist_mat)
  if (is.null(rn)) {
    rownames(hist_mat) <- as.character(meta_df$site)
    rn <- rownames(hist_mat)
  }
  idx <- match(rn, as.character(meta_df$site))
  if (any(is.na(idx))) abort("occu_build_umf(): rownames(hist_mat) must match meta_df$site.")
  meta_use <- meta_df[idx, , drop = FALSE]
  if (!all(rownames(hist_mat) == as.character(meta_use$site))) abort("occu_build_umf(): internal alignment failure.")
  type <- factor(as.character(meta_use$type), levels = c("dry","wet"))
  unmarked::unmarkedFrameOccu(y = hist_mat, siteCovs = data.frame(type = type))
}

# 10) fit simple occupancy model ψ(~ type), p(~1)
fit_simple_occu <- function(umf, with_type = TRUE) {
  if (!requireNamespace("unmarked", quietly = TRUE)) {
    abort("Package 'unmarked' is required for the occupancy demo. Add it via dependencies.R then renv::restore().")
  }
  warn_msgs <- character(0)
  res <- withCallingHandlers({
    # Model: ψ ~ type (if with_type) else ψ ~ 1; p ~ 1
    m <- if (isTRUE(with_type)) unmarked::occu(~ 1 ~ type, data = umf) else unmarked::occu(~ 1 ~ 1, data = umf)

    if (isTRUE(with_type)) {
      # Predict ψ by habitat
      newd <- data.frame(type = factor(c("dry","wet"), levels = c("dry","wet")))
      ps  <- unmarked::predict(m, type = "state", newdata = newd)
      tidy <- tibble(
        param = c("psi_dry","psi_wet"),
        est   = as.numeric(ps$Predicted),
        lwr   = as.numeric(ps$lower),
        upr   = as.numeric(ps$upper),
        link  = "logit^-1",
        method = "Wald"
      )
    } else {
      # ψ overall via backTransform if available; otherwise predict()
      bt_ok <- try(unmarked::backTransform(m, type = "state"), silent = TRUE)
      if (!inherits(bt_ok, "try-error")) {
        ci <- try(stats::confint(bt_ok), silent = TRUE)
        est <- as.numeric(bt_ok@estimate)
        if (!inherits(ci, "try-error") && length(ci) >= 2) {
          lwr <- as.numeric(ci[1]); upr <- as.numeric(ci[2])
        } else {
          pr <- unmarked::predict(m, type = "state")
          est <- as.numeric(pr$Predicted[1]); lwr <- as.numeric(pr$lower[1]); upr <- as.numeric(pr$upper[1])
        }
      } else {
        pr <- unmarked::predict(m, type = "state")
        est <- as.numeric(pr$Predicted[1]); lwr <- as.numeric(pr$lower[1]); upr <- as.numeric(pr$upper[1])
      }
      tidy <- tibble(param = "psi_overall", est = est, lwr = lwr, upr = upr, link = "logit^-1", method = "Wald")
    }

    # Detection probability p (intercept-only)
    pd <- unmarked::predict(m, type = "det")
    tidy <- bind_rows(tidy, tibble(param = "p", est = as.numeric(pd$Predicted[1]), lwr = as.numeric(pd$lower[1]), upr = as.numeric(pd$upper[1]), link = "logit^-1", method = "Wald"))
    # Scale labelling to prevent confusion
    tidy$scale <- ifelse(grepl("^psi", tidy$param), "psi (site-level)", "p (per-day)")

    # Figure: ψ points
    if (isTRUE(with_type)) {
      dfp <- tidy[tidy$param %in% c("psi_dry","psi_wet"), , drop = FALSE]
      dfp$habitat <- c("dry","wet")
      fig <- ggplot(dfp, aes(x = habitat, y = est)) +
        ggplot2::geom_pointrange(aes(ymin = lwr, ymax = upr)) +
        theme_minimal() +
        labs(title = "Occupancy ψ by habitat (95% CI)", x = NULL, y = "ψ")
    } else {
      dfp <- tidy[tidy$param == "psi_overall", , drop = FALSE]
      dfp$label <- "overall"
      fig <- ggplot(dfp, aes(x = label, y = est)) +
        ggplot2::geom_pointrange(aes(ymin = lwr, ymax = upr)) +
        theme_minimal() +
        labs(title = "Occupancy ψ (overall, 95% CI)", x = NULL, y = "ψ")
    }

    list(model = m, tidy = tidy, fig = fig)
  }, warning = function(w) {
    warn_msgs <<- c(warn_msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })

  note <- if (length(warn_msgs)) paste(unique(warn_msgs), collapse = " | ") else ""
res$tidy <- dplyr::mutate(res$tidy, convergence_note = note)
  if (nzchar(note)) res$fig <- res$fig + ggplot2::labs(caption = note)
  res
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# 11) long-form detection history for inspection (site × day rows)
detection_history_long <- function(hist_mat, meta_df) {
  if (!is.matrix(hist_mat)) abort("detection_history_long(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("detection_history_long(): meta_df must have site and type.")
  df <- melt(hist_mat, varnames = c("site","day"), value.name = "det")
  df$day <- as.integer(sub("^d", "", as.character(df$day)))
  meta_small <- distinct(meta_df[, c("site","type")])
  left_join(df, meta_small, by = "site")
}

# 12) per-site summary for occupancy inspection
summarise_occ_site <- function(hist_mat, meta_df) {
  if (!is.matrix(hist_mat)) abort("summarise_occ_site(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("summarise_occ_site(): meta_df must have site and type.")
  rn <- rownames(hist_mat)
  if (is.null(rn)) rn <- meta_df$site
  det_days <- rowSums(hist_mat > 0, na.rm = TRUE)
  first_day <- apply(hist_mat, 1, function(v) {
    i <- which(v > 0)
    if (length(i)) min(i) else NA_integer_
  })
  last_day <- apply(hist_mat, 1, function(v) {
    i <- which(v > 0)
    if (length(i)) max(i) else NA_integer_
  })
  hist_str <- apply(hist_mat, 1, function(v) paste(as.integer(v), collapse = ""))
  out <- tibble(
    site = rn,
    n_days = ncol(hist_mat),
    n_detect_days = as.integer(det_days),
    ever_detected = det_days > 0,
    first_day = as.integer(first_day),
    last_day  = as.integer(last_day),
    history   = hist_str
  )
  left_join(out, distinct(meta_df[, c("site","type")]), by = "site")
}
# 11) detection-rate summary for picking a demo species
detection_rate_summary <- function(hist_mat, meta_df) {
  if (!is.matrix(hist_mat)) abort("detection_rate_summary(): hist_mat must be a matrix.")
  if (!all(c("site","type") %in% names(meta_df))) abort("detection_rate_summary(): meta_df must have site and type.")
  sp <- attr(hist_mat, "species") %||% "(unknown)"
  ever <- rowSums(hist_mat > 0) > 0
  overall <- mean(ever)
  # by habitat
  rn <- rownames(hist_mat)
  if (is.null(rn)) rn <- meta_df$site
  meta_small <- distinct(meta_df[, c("site","type")])
  meta_use <- meta_small[match(rn, meta_small$site), , drop = FALSE]
  by_type <- split(ever, as.character(meta_use$type))
  prop_dry <- if (!is.null(by_type$dry)) mean(by_type$dry) else NA_real_
  prop_wet <- if (!is.null(by_type$wet)) mean(by_type$wet) else NA_real_
  tibble(
    species = sp,
    n_sites = length(ever),
    n_detect_sites = sum(ever),
    prop_overall = overall,
    prop_dry = prop_dry,
    prop_wet = prop_wet
  )
}

# 12) choose a pedagogical species (closest to 0.5 overall, prefer 0.2–0.8)
choose_occ_species <- function(summary_tbl) {
  if (!all(c("species","prop_overall") %in% names(summary_tbl))) abort("choose_occ_species(): summary_tbl missing columns.")
  df <- summary_tbl
  df$score <- abs(df$prop_overall - 0.5)
  # prefer those within 0.2–0.8; else take closest overall
  in_band <- df$prop_overall >= 0.2 & df$prop_overall <= 0.8
  cand <- if (any(in_band)) df[in_band, , drop = FALSE] else df
  cand <- cand[order(cand$score, decreasing = FALSE), , drop = FALSE]
  as.character(cand$species[1])
}
