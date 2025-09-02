# R/helpers_interactions.R — co-occurrence diagnostics under independence

import::here(abort, .from = "rlang")
import::here(tibble, .from = "tibble")
import::here(dplyr_select = select, left_join, group_by, summarise, n, ungroup, mutate, arrange, distinct, .from = "dplyr")
import::here(pivot_longer, pivot_wider, .from = "tidyr")
import::here(melt, .from = "reshape2")

# 1) presence–absence from community matrix
pa_from_comm <- function(comm_mat) {
  mat <- as.matrix(comm_mat)
  if (is.null(rownames(mat))) rownames(mat) <- paste0("S", seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("sp", seq_len(ncol(mat)))
  pa <- ifelse(mat > 0, 1L, 0L)
  storage.mode(pa) <- "integer"
  pa
}

# 2) species selection for co-occurrence
select_coocc_species <- function(pa_mat, include = character(), k_max = 16L, min_sites = 4L) {
  mat <- as.matrix(pa_mat)
  occ <- colSums(mat, na.rm = TRUE)
  # pre-include always kept
  inc <- intersect(include, colnames(mat))
  # rank by occupancy
  ranked <- names(sort(occ, decreasing = TRUE))
  chosen <- unique(c(inc, ranked))
  # filter by min_sites and cap
  ok <- chosen[occ[chosen] >= as.integer(min_sites)]
  out <- head(ok, as.integer(k_max))
  if (length(out) < 2L) abort("select_coocc_species(): fewer than 2 species satisfy min_sites.")
  out
}

# 3) pairwise stats under fixed-margins independence (optionally by group)
coocc_pair_stats <- function(pa_mat, meta_df = NULL, group_col = NULL) {
  mat <- as.matrix(pa_mat)
  if (!is.null(group_col)) {
    if (is.null(meta_df) || !("site" %in% names(meta_df)) || !(group_col %in% names(meta_df))) {
      abort("coocc_pair_stats(): meta_df must contain site and group_col when group_col is provided.")
    }
    labs <- rownames(mat)
    if (is.null(labs)) labs <- meta_df$site
    # align meta to rows
    meta_use <- meta_df[match(labs, meta_df$site), , drop = FALSE]
    if (anyNA(meta_use$site)) abort("coocc_pair_stats(): cannot align meta to PA matrix rows.")
    groups <- unique(as.character(meta_use[[group_col]]))
    out <- lapply(groups, function(g) {
      idx <- which(as.character(meta_use[[group_col]]) == g)
      if (length(idx) < 3L) return(NULL)
      sub <- mat[idx, , drop = FALSE]
      tib <- coocc_pair_stats(sub)
      if (nrow(tib)) tib$group <- g
      tib
    })
    out <- Filter(Negate(is.null), out)
    if (!length(out)) return(tibble())
    res <- do.call(rbind, out)
    # ensure group first col
    res <- res[, c("group", setdiff(names(res), "group"))]
    return(res)
  }

  sp <- colnames(mat)
  N <- nrow(mat)
  if (N < 3L) abort("coocc_pair_stats(): need at least 3 sites.")
  Xi <- colSums(mat, na.rm = TRUE)
  pairs <- combn(seq_along(sp), 2L)
  res <- apply(pairs, 2L, function(ij) {
    i <- ij[1]; j <- ij[2]
    xi <- Xi[i]; xj <- Xi[j]
    o  <- sum(mat[, i] & mat[, j])
    mu <- xi * xj / N
    var <- (xi * xj * (N - xi) * (N - xj)) / (N^2 * (N - 1))
    z <- if (var > 0) (o - mu) / sqrt(var) else 0
    p <- 2 * stats::pnorm(-abs(z))
    tibble(
      sp_i = sp[i], sp_j = sp[j], N = as.integer(N), Xi = as.integer(xi), Xj = as.integer(xj),
      O = as.integer(o), E = as.numeric(mu), z = as.numeric(z), p_norm = as.numeric(p)
    )
  })
  df <- do.call(rbind, res)
  # BH-FDR q-values as a reality check (no plotting change)
  df$q_bh <- stats::p.adjust(df$p_norm, method = "BH")
  df
}

# 4) build symmetric z-matrix (diag = NA)
coocc_z_matrix <- function(stats_tbl, group = "all") {
  df <- stats_tbl
  if (!nrow(df)) abort("coocc_z_matrix(): empty stats table.")
  if ("group" %in% names(df)) df <- df[df$group == group, , drop = FALSE]
  spp <- unique(c(as.character(df$sp_i), as.character(df$sp_j)))
  zmat <- matrix(NA_real_, nrow = length(spp), ncol = length(spp), dimnames = list(spp, spp))
  for (k in seq_len(nrow(df))) {
    i <- as.character(df$sp_i[k]); j <- as.character(df$sp_j[k]); z <- as.numeric(df$z[k])
    zmat[i, j] <- z
    zmat[j, i] <- z
  }
  diag(zmat) <- NA_real_
  zmat
}

# helper: z-matrices per group
coocc_z_by_group <- function(stats_tbl) {
  if (!("group" %in% names(stats_tbl))) abort("coocc_z_by_group(): 'group' column expected.")
  grps <- unique(as.character(stats_tbl$group))
  out <- lapply(grps, function(g) coocc_z_matrix(stats_tbl, group = g))
  names(out) <- grps
  out
}
