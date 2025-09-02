# R/helpers_activity.R — diel activity overlap (Δ̂) and visuals

import::here(abort, .from = "rlang")
import::here(tibble, .from = "tibble")
import::here(filter, mutate, left_join, group_by, summarise, .from = "dplyr")
import::from("lubridate", with_tz, hour, minute, second)
import::from("overlap", overlapEst, bootstrap)
import::from("circular", circular, density.circular)
import::here(ggplot, aes, geom_line, geom_ribbon, facet_wrap, theme_minimal, labs, scale_y_continuous, scale_x_continuous, .from = "ggplot2")

.rad_from_datetime <- function(dt) {
  # dt assumed POSIXct in any tz; compute Hobart-local time-of-day in radians
  loc <- with_tz(dt, tzone = "Australia/Hobart")
  s <- as.numeric(hour(loc)) * 3600 + as.numeric(minute(loc)) * 60 + as.numeric(second(loc))
  (s / 86400) * 2 * pi
}

# 6) times by species (optionally by habitat type)
times_by_species <- function(events_df, species, type = NULL) {
  cols <- c("common","datetime")
  if (!all(cols %in% names(events_df))) abort("times_by_species(): events_df must have 'common' and 'datetime'.")
  df <- events_df[events_df$common == species, , drop = FALSE]
  if (!is.null(type) && ("type" %in% names(events_df))) {
    df <- df[df$type == type, , drop = FALSE]
  }
  if (nrow(df) == 0) return(numeric(0))
  .rad_from_datetime(df$datetime)
}

# 7) estimator choice rule
choose_overlap_estimator <- function(n1, n2) {
  if (min(as.integer(n1), as.integer(n2)) >= 75L) "Dhat4" else "Dhat1"
}

# 8) overlap Δ̂ and bootstrap CI by pair and habitat groups
overlap_by_pair <- function(events_df, pair, groups = c("wet","dry"), n_boot = 500L, seed = 1L) {
  # robustify pair handling for list/character vectors
  if (is.list(pair) && length(pair) == 1L) pair <- pair[[1]]
  if (is.list(pair) && length(pair) == 2L && all(vapply(pair, length, 1L) == 1L)) pair <- c(pair[[1]], pair[[2]])
  if (length(pair) != 2L) abort("overlap_by_pair(): pair must have length 2.")
  sp1 <- as.character(pair[[1]]); sp2 <- as.character(pair[[2]])

  out <- lapply(groups, function(g) {
    t1 <- times_by_species(events_df, sp1, type = g)
    t2 <- times_by_species(events_df, sp2, type = g)
    n1 <- length(t1); n2 <- length(t2)
    est <- choose_overlap_estimator(n1, n2)
    # compute Δ̂ and bootstrap CI
    if (n1 == 0 || n2 == 0) {
      delta <- NA_real_; ci <- c(NA_real_, NA_real_)
      p_watson <- NA_real_
    } else {
      delta <- overlap::overlapEst(t1, t2, type = est)
      set.seed(as.integer(seed))
      boots <- overlap::bootstrap(t1, t2, as.integer(n_boot), type = est)
      ci <- stats::quantile(boots, c(0.025, 0.975), na.rm = TRUE)
      # Watson two-sample test (if available)
      p_watson <- NA_real_
      if ("watson.two.test" %in% getNamespaceExports("circular")) {
        c1 <- circular::circular(t1, units = "radians", template = "clock24")
        c2 <- circular::circular(t2, units = "radians", template = "clock24")
        wt <- try(circular::watson.two.test(c1, c2), silent = TRUE)
        if (!inherits(wt, "try-error") && !is.null(wt$p.value)) p_watson <- as.numeric(wt$p.value)
      }
    }
    tibble(
      species1 = sp1, species2 = sp2, type = g, n1 = as.integer(n1), n2 = as.integer(n2),
      estimator = est, delta_hat = as.numeric(delta), lwr = as.numeric(ci[[1]]), upr = as.numeric(ci[[2]]),
      p_watson = p_watson,
      n_boot = as.integer(n_boot), seed = as.integer(seed)
    )
  })
  do.call(rbind, out)
}

# 9) circular density to 24h grid
density_circular_df <- function(times, grid = 256L) {
  if (length(times) == 0) return(tibble(hour = numeric(0), density = numeric(0)))
  if (length(times) < 2) {
    # fallback: near-constant low density to avoid bw failure
    hour <- seq(0, 24, length.out = as.integer(grid))
    return(tibble(hour = hour, density = rep(1 / 24, length(hour))))
  }
  circ <- suppressWarnings(circular::circular(times, units = "radians", template = "clock24", modulo = "2pi", zero = 0, rotation = "counter"))
  den <- try(suppressWarnings(circular::density.circular(circ, bw = NULL, kernel = "vonmises", from = 0, to = 2 * pi, n = as.integer(grid))), silent = TRUE)
  if (inherits(den, "try-error")) {
    den <- try(suppressWarnings(circular::density.circular(circ, bw = 10, kernel = "vonmises", from = 0, to = 2 * pi, n = as.integer(grid))), silent = TRUE)
  }
  if (inherits(den, "try-error")) {
    hour <- seq(0, 24, length.out = as.integer(grid))
    return(tibble(hour = hour, density = rep(1 / 24, length(hour))))
  }
  hour <- (as.numeric(den$x) / (2 * pi)) * 24
  tibble(hour = hour, density = as.numeric(den$y))
}

# 10) plot activity overlap with shaded common area
plot_activity_overlap <- function(events_df, pair, groups = c("wet","dry")) {
  # flatten pair input like in overlap_by_pair
  if (is.list(pair) && length(pair) == 1L) pair <- pair[[1]]
  if (is.list(pair) && length(pair) == 2L && all(vapply(pair, length, 1L) == 1L)) pair <- c(pair[[1]], pair[[2]])
  if (length(pair) != 2L) abort("plot_activity_overlap(): pair must have length 2.")
  sp1 <- as.character(pair[[1]]); sp2 <- as.character(pair[[2]])
  # build per-group densities
  df_list <- lapply(groups, function(g) {
    t1 <- times_by_species(events_df, sp1, type = g)
    t2 <- times_by_species(events_df, sp2, type = g)
    d1 <- density_circular_df(t1)
    d2 <- density_circular_df(t2)
    if (nrow(d1) == 0 || nrow(d2) == 0) {
      dens <- tibble(hour = numeric(0), density1 = numeric(0), density2 = numeric(0), type = character(0))
    } else {
      if (nrow(d1) == nrow(d2) && max(abs(d1$hour - d2$hour)) < 1e-8) {
        dens <- tibble(hour = d1$hour, density1 = d1$density, density2 = d2$density, type = g)
      } else {
        # de-duplicate hours then join with explicit many-to-many relationship (expected on rounded floats)
        d1u <- dplyr::distinct(d1, hour, .keep_all = TRUE)
        d2u <- dplyr::distinct(d2, hour, .keep_all = TRUE)
        tmp <- dplyr::left_join(d1u, d2u, by = "hour", relationship = "many-to-many", suffix = c("1","2"))
        colnames(tmp) <- c("hour","density1","density2")
        tmp$type <- g
        dens <- tmp
      }
    }
    dens
  })
  dens_all <- do.call(rbind, df_list)

  # summary for titles
  summ <- overlap_by_pair(events_df, pair, groups = groups, n_boot = 200L, seed = 1L)
  summ$title <- sprintf("Δ̂ = %.2f (95%% CI %.2f–%.2f), n=(%d,%d), est=%s",
                        summ$delta_hat, summ$lwr, summ$upr, summ$n1, summ$n2, summ$estimator)
  low_n <- summ$n1 < 10 | summ$n2 < 10
  summ$subtitle <- ifelse(low_n, "low n—interpret cautiously", "")

  # build plot
  if (nrow(dens_all) == 0) {
    return(ggplot2::ggplot() + theme_minimal() + labs(title = paste(sp1, "vs", sp2), subtitle = "no data"))
  }
  dens_all$overlap <- pmin(dens_all$density1, dens_all$density2)
  p <- ggplot(dens_all, aes(x = hour)) +
    geom_ribbon(aes(ymin = 0, ymax = overlap), fill = "grey80") +
    geom_line(aes(y = density1, colour = pair[1])) +
    geom_line(aes(y = density2, colour = pair[2])) +
    facet_wrap(~ type, ncol = 1, scales = "free_y") +
    theme_minimal() +
    labs(title = paste(sp1, "vs", sp2, "activity overlap"), x = "Hour (local)", y = "Density")

  # add facet titles with Δ̂
  # Note: direct per-facet titles require strip labels; we keep overall title and list Δ̂ in caption
  cap <- paste(sprintf("%s: %s", summ$type, summ$title), collapse = " | ")
  if (any(nzchar(summ$subtitle))) cap <- paste(cap, "; notes: low sample size", sep = " ")
  p + labs(caption = cap)
}
