# R/helpers_beta.R — beta diversity distances and inference

import::here(vegdist, adonis2, betadisper, permutest, metaMDS, scores, .from = "vegan")
import::here(tibble, distinct, left_join, mutate, .from = "dplyr")
import::here(abort, .from = "rlang")
import::here(local_seed, .from = "withr")

#' Compute a between-site distance object (Bray–Curtis or Jaccard)
#'
#' @param comm_mat site × species abundance matrix (event counts)
#' @param method   one of "bray" (abundance) or "jaccard" (presence–absence)
#' @return stats::dist object with site labels
#' @export
compute_dist <- function(comm_mat, method = c("bray", "jaccard")) {
  method <- match.arg(method)

  mat <- as.matrix(comm_mat)
  if (nrow(mat) < 2L) abort("compute_dist(): need >= 2 sites (rows) to compute distances.")

  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("S", seq_len(nrow(mat)))
  }

  if (method == "bray") {
    d <- vegdist(mat, method = "bray")
  } else {
    # presence–absence for Jaccard
    d <- vegdist(mat, method = "jaccard", binary = TRUE)
  }
  d
}

#' PERMANOVA + dispersion test (betadisper) for a single grouping factor
#'
#' @param dist_obj  a stats::dist object
#' @param meta_df   data.frame with columns site and the grouping column
#' @param group_col string name of grouping column in meta_df (default "type")
#' @param n_perm    number of permutations for tests (default 999)
#' @param seed      RNG seed for permutation reproducibility
#' @return tibble(distance, term, df, pseudo_F, R2, p_perm, disp_F, disp_p, note)
#' @export
permanova_and_dispersion <- function(dist_obj, meta_df, group_col = "type",
                                     n_perm = 999L, seed = 1L) {
  if (!("site" %in% names(meta_df))) abort("permanova_and_dispersion(): meta_df must contain 'site'.")
  if (!(group_col %in% names(meta_df))) abort(sprintf("permanova_and_dispersion(): meta_df must contain '%s'.", group_col))

  # align metadata to distance labels
  labs <- attr(dist_obj, "Labels")
  if (is.null(labs)) labs <- attr(dist_obj, "labels")
  if (is.null(labs)) labs <- rownames(as.matrix(dist_obj))
  if (is.null(labs)) abort("permanova_and_dispersion(): cannot retrieve site labels from dist object.")

  meta_small <- distinct(meta_df, site, !!group_col)
  # handle dynamic column without tidy-eval dependency: use base indexing
  meta_small <- meta_df[, c("site", group_col)]
  meta_small <- unique(meta_small)
  colnames(meta_small) <- c("site", "group")

  meta_use <- meta_small[match(labs, meta_small$site), , drop = FALSE]
  if (anyNA(meta_use$group)) abort("permanova_and_dispersion(): some sites in the distance object are missing from metadata.")

  # determine method label from dist attributes
  method <- attr(dist_obj, "method")
  if (is.null(method)) method <- "unknown"

  local_seed(as.integer(seed))

  # PERMANOVA (adonis2)
  ad <- adonis2(dist_obj ~ group, data = meta_use, permutations = as.integer(n_perm))
  # adonis2 returns a data.frame-like; row name for the term is 'group'
  if (!("R2" %in% names(ad))) abort("permanova_and_dispersion(): adonis2 output missing R2.")
  if (!("F" %in% names(ad))) abort("permanova_and_dispersion(): adonis2 output missing F.")

  # row for the grouping term is usually first; find by rowname if present
  rn <- rownames(ad)
  idx <- if (!is.null(rn)) which(rn == "group") else 1L
  if (length(idx) == 0L) idx <- 1L
  term_df <- ad[idx, , drop = FALSE]

  R2 <- as.numeric(term_df$R2)
  Fv <- as.numeric(term_df$F)
  p  <- suppressWarnings(as.numeric(term_df$`Pr(>F)`))
  Df <- suppressWarnings(as.integer(term_df$Df))

  if (!is.na(R2) && (R2 < 0 || R2 > 1)) abort("permanova_and_dispersion(): R2 outside [0,1].")
  if (!is.na(p)  && (p  < 0 || p  > 1)) abort("permanova_and_dispersion(): p-value outside [0,1].")

  # Dispersion: requires >= 2 sites per group
  counts <- table(meta_use$group)
  small <- any(counts < 2L)
  if (small) {
    disp_F <- NA_real_
    disp_p <- NA_real_
    note   <- "Group size <2"
  } else {
    bd  <- betadisper(dist_obj, meta_use$group)
    pt  <- permutest(bd, permutations = as.integer(n_perm))
    disp_F <- as.numeric(pt$tab["Groups", "F"])
    disp_p <- as.numeric(pt$tab["Groups", "Pr(>F)"])
    note <- if (!is.na(disp_p) && disp_p < 0.05) "Dispersion differs (interpret location cautiously)" else ""
  }

  tibble(
    distance = method,
    term     = group_col,
    df       = Df,
    pseudo_F = Fv,
    R2       = R2,
    p_perm   = p,
    disp_F   = disp_F,
    disp_p   = disp_p,
    note     = note
  )
}

