# R/helpers_glm.R — effort-offset GLMs for habitat selection

import::here(select, mutate, transmute, rename, left_join, inner_join, distinct,
             filter, group_by, summarise, ungroup, across, bind_cols, tibble,
             .from = "dplyr")
import::here(pivot_wider, .from = "tidyr")
import::here(abort, .from = "rlang")
import::from("stats", "glm", "poisson", "quasipoisson", "model.matrix", "coef",
             "vcov", "predict", "qnorm", "AIC")
import::here(ggplot, aes, geom_point, geom_errorbar, geom_line, geom_ribbon,
             theme_minimal, labs, scale_x_continuous, scale_y_continuous,
             scale_x_log10, geom_vline, .from = "ggplot2")

z_score <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

#' Prepare per-species GLM data with offsets and scaled covariates
#'
#' @param events_df events table: site, common, event
#' @param site_df   site metadata: site, type, covariates
#' @param effort_df effort table: site, trap_nights
#' @param species   character scalar
#' @param covars    character vector including "type" and numeric covariates
#' @return tibble with columns: site, species, y, trap_nights, type, *_z; attr("dropped_sites") recorded
#' @export
prepare_species_glm_data <- function(events_df, site_df, effort_df, species,
                                     covars = c("type","leaf_litter","cwd","can_foliage")) {
  req_ev <- c("site","common","event")
  if (!all(req_ev %in% names(events_df))) abort("prepare_species_glm_data(): events_df missing required columns.")
  if (!all(c("site","type") %in% names(site_df))) abort("prepare_species_glm_data(): site_df missing 'site' or 'type'.")
  if (!all(c("site","trap_nights") %in% names(effort_df))) abort("prepare_species_glm_data(): effort_df missing 'site' or 'trap_nights'.")
  if (!("type" %in% covars)) abort("prepare_species_glm_data(): covars must include 'type'.")

  # species counts per site
  y_tbl <- events_df |>
    filter(common == species) |>
    group_by(site, common, event) |>
    summarise(n = 1L, .groups = "drop") |>
    group_by(site) |>
    summarise(y = n(), .groups = "drop")

  # all sites (including zero counts)
  all_sites <- distinct(site_df, site)
  y_all <- left_join(all_sites, y_tbl, by = "site")
  y_all$y[is.na(y_all$y)] <- 0L

  # bring effort and covariates
  meta <- left_join(site_df, select(effort_df, site, trap_nights), by = "site")

  # filter to available covariates
  keep_cols <- unique(c("site","type","trap_nights", setdiff(covars, "type")))
  meta <- meta[, intersect(keep_cols, names(meta)), drop = FALSE]

  dat <- left_join(y_all, meta, by = "site")

  # ensure factor baseline for type: dry then wet
  if ("type" %in% names(dat)) dat$type <- factor(dat$type, levels = c("dry","wet"))

  # z-score numeric covariates
  num_covars <- setdiff(covars, "type")
  for (v in num_covars) {
    zname <- paste0(v, "_z")
    dat[[zname]] <- z_score(dat[[v]])
  }

  # drop rows with missing covariates or invalid effort
  z_cols <- paste0(setdiff(covars, "type"), "_z")
  to_drop <- !is.finite(dat$trap_nights) | dat$trap_nights <= 0
  for (v in z_cols) to_drop <- to_drop | !is.finite(dat[[v]])
  if ("type" %in% covars) to_drop <- to_drop | is.na(dat$type)

  dropped <- sum(to_drop, na.rm = TRUE)
  dat <- dat[!to_drop, , drop = FALSE]

  if (nrow(dat) < 4L) abort(sprintf("prepare_species_glm_data(): < 4 sites remain for %s.", species))

  dat$species <- species

  attr(dat, "dropped_sites") <- dropped
  attr(dat, "species") <- species
  dat
}

#' Fit effort-offset GLM with Poisson/quasi-Poisson/NB selection
#'
#' @param df prepared data from prepare_species_glm_data()
#' @param species species name (character)
#' @return list(model, family, dispersion, aic, theta, n_sites, n_zero_sites, dropped_sites, note)
#' @export
fit_species_glm <- function(df, species) {
  if (!all(c("y","trap_nights") %in% names(df))) abort("fit_species_glm(): df missing y or trap_nights.")
  # decide terms
  has_type <- "type" %in% names(df) && length(unique(stats::na.omit(df$type))) >= 2L
  note_vec <- character(0)
  # Complete separation guard: all events in a single level of type
  if (has_type) {
    sums_by_type <- tapply(df$y, df$type, sum, na.rm = TRUE)
    if (length(sums_by_type) >= 2L && any(sums_by_type == 0)) {
      has_type <- FALSE
      note_vec <- c(note_vec, "type dropped (events only in one level; separation)")
    }
  }
  terms <- c(if (has_type) "type" else NULL, "leaf_litter_z", "cwd_z", "can_foliage_z")
  terms <- terms[terms %in% names(df)]
  fml <- stats::as.formula(paste("y ~", paste(terms, collapse = " + ")))

  # Poisson base fit
  m_p <- glm(fml, family = poisson(link = "log"), offset = log(trap_nights), data = df, model = TRUE)
  dfres <- max(1, m_p$df.residual)
  phi <- m_p$deviance / dfres

  avg_traps <- mean(df$trap_nights, na.rm = TRUE)
  if (!is.finite(avg_traps) || is.na(avg_traps)) avg_traps <- 1

  chosen <- list(model = m_p, family = "poisson", dispersion = phi,
                  aic = AIC(m_p), theta = NA_real_,
                  n_sites = nrow(df), n_zero_sites = sum(df$y == 0),
                 dropped_sites = attr(df, "dropped_sites") %||% 0L,
                 avg_traps = avg_traps,
                 species = as.character(species),
                 note = {
                   msg <- if (!has_type) c("type dropped (single level)") else character(0)
                   paste(c(note_vec, msg), collapse = "; ")
                 })

  if (phi <= 1.5) return(chosen)

  # Try NB
  m_nb <- try(MASS::glm.nb(fml, offset = log(trap_nights), data = df, model = TRUE), silent = TRUE)
  if (!inherits(m_nb, "try-error")) {
    aic_nb <- AIC(m_nb)
    if (is.finite(aic_nb) && (aic_nb + 2 < chosen$aic)) {
      phi_nb <- m_nb$deviance / max(1, m_nb$df.residual)
      return(list(model = m_nb, family = "negbin", dispersion = phi_nb, aic = aic_nb,
                  theta = as.numeric(m_nb$theta), n_sites = chosen$n_sites,
                  n_zero_sites = chosen$n_zero_sites, dropped_sites = chosen$dropped_sites,
                  avg_traps = avg_traps,
                  species = as.character(species),
                  note = chosen$note))
    }
  }

  # Fallback: quasi-Poisson (robust SEs)
  m_q <- glm(fml, family = quasipoisson(link = "log"), offset = log(trap_nights), data = df, model = TRUE)
  phi_q <- m_q$deviance / max(1, m_q$df.residual)
  list(model = m_q, family = "quasipoisson", dispersion = phi_q, aic = NA_real_, theta = NA_real_,
       n_sites = chosen$n_sites, n_zero_sites = chosen$n_zero_sites,
       dropped_sites = chosen$dropped_sites, avg_traps = avg_traps,
       species = as.character(species), note = chosen$note)
}

#' IRRs with 95% Wald CIs from a fitted model info list
#'
#' @param mod_info list from fit_species_glm()
#' @return tibble(species, term, irr, lwr, upr, p_value, family)
#' @export
tidy_irrs <- function(mod_info) {
  mod <- mod_info$model
  cf  <- summary(mod)$coefficients
  terms <- rownames(cf)
  est <- cf[, 1]
  se  <- cf[, 2]
  z   <- 1.96
  irr <- exp(est)
  lwr <- exp(est - z * se)
  upr <- exp(est + z * se)
  pv  <- if (ncol(cf) >= 4) cf[, 4] else NA_real_
  tibble(
    species = as.character(mod_info$species %||% NA_character_),
    term = terms,
    irr = as.numeric(irr),
    lwr = as.numeric(lwr),
    upr = as.numeric(upr),
    p_value = as.numeric(pv),
    family = mod_info$family
  )
}

#' Diagnostics table
#'
#' @param mod_info list from fit_species_glm()
#' @return tibble(species, family, dispersion, aic, theta, n_sites, n_zero_sites, dropped_sites, note)
#' @export
diagnostics_table <- function(mod_info) {
  tibble(
    species = as.character(mod_info$species %||% NA_character_),
    family = mod_info$family,
    dispersion = as.numeric(mod_info$dispersion),
    n_sites = as.integer(mod_info$n_sites),
    n_zero_sites = as.integer(mod_info$n_zero_sites),
    dropped_sites = as.integer(mod_info$dropped_sites),
    note = mod_info$note %||% ""
  )
}

mm_and_eta <- function(mod, newdata) {
  tt <- stats::terms(mod)
  tt_x <- stats::delete.response(tt)
  mm <- stats::model.matrix(tt_x, data = newdata,
                            contrasts.arg = mod$contrasts, xlev = mod$xlevels,
                            na.action = stats::na.pass)
  eta <- as.vector(mm %*% stats::coef(mod))
  # Always add offset using trap_nights from newdata if present
  off <- if ("trap_nights" %in% names(newdata)) log(newdata$trap_nights) else rep(0, nrow(newdata))
  eta <- eta + off
  list(mm = mm, eta = eta)
}

#' Partial dependence (Wald bands) keeping other covariates at 0 (their z-means)
#'
#' @param mod_info list from fit_species_glm()
#' @param df optional original df (falls back to model frame)
#' @param var variable name: "type" or *_z
#' @param grid_n number of grid points for numeric variables
#' @return ggplot
#' @export
partial_dependence <- function(mod_info, df = NULL, var, grid_n = 50L) {
  mod <- mod_info$model
  if (is.null(mod) || !inherits(mod, "glm")) {
    abort("partial_dependence(): mod_info$model is not a fitted 'glm'.")
  }
  if (is.null(df)) df <- mod$model
  if (is.null(df)) abort("partial_dependence(): cannot locate training data from model.")

  # Baseline settings derived from the model frame to match types/levels
  avg_traps <- mod_info$avg_traps %||% 1
  base <- data.frame(trap_nights = avg_traps)
  mf <- try(stats::model.frame(mod), silent = TRUE)
  if (!inherits(mf, "try-error")) {
    resp_idx <- attr(stats::terms(mod), "response")
    preds <- colnames(mf)
    if (!is.null(resp_idx) && resp_idx > 0 && resp_idx <= length(preds)) {
      preds <- preds[-resp_idx]
    }
    for (nm in preds) {
      if (is.factor(mf[[nm]])) {
        lv <- levels(mf[[nm]])
        base[[nm]] <- factor(if ("dry" %in% lv) "dry" else lv[1], levels = lv)
      } else {
        base[[nm]] <- 0
      }
    }
  }
  # Ensure any factor terms represented in coefficients are present (e.g., type)
  coef_terms <- names(stats::coef(mod))
  if (any(grepl("^type", coef_terms)) && !("type" %in% names(base))) {
    lv <- if (!is.null(mod$xlevels) && "type" %in% names(mod$xlevels)) mod$xlevels$type else c("dry","wet")
    base$type <- factor(if ("dry" %in% lv) "dry" else lv[1], levels = lv)
  }
  if (!("leaf_litter_z" %in% names(base))) base$leaf_litter_z <- 0
  if (!("cwd_z" %in% names(base))) base$cwd_z <- 0
  if (!("can_foliage_z" %in% names(base))) base$can_foliage_z <- 0

  if (identical(var, "type")) {
    # Two points: dry vs wet; include 'type' even if not used in the model
    nd <- base[rep(1, times = 2L), , drop = FALSE]
    if (!inherits(mf, "try-error") && "type" %in% names(mf)) {
      lv <- levels(mf$type)
      if (length(lv) >= 2L) nd$type <- factor(lv[1:2], levels = lv) else nd$type <- factor(rep(lv[1], 2L), levels = lv)
    } else {
      lv <- if (!is.null(mod$xlevels) && "type" %in% names(mod$xlevels)) mod$xlevels$type else c("dry","wet")
      sel <- if (length(lv) >= 2L) lv[1:2] else rep(lv[1], 2L)
      nd$type <- factor(sel, levels = lv)
    }
    pred <- predict_wald(mod, nd)
    nd$fit <- pred$fit; nd$lwr <- pred$lwr; nd$upr <- pred$upr
    ggplot(nd, aes(x = type, y = fit)) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.15) +
      theme_minimal() +
      labs(title = as.character(mod_info$species), y = "Predicted events (offset = mean trap-nights)", x = "Habitat type")
  } else {
    # numeric var in z-scale; be robust to missing/constant columns
    vals <- if (var %in% names(df)) suppressWarnings(as.numeric(df[[var]])) else numeric(0)
    if (length(vals) == 0L || all(!is.finite(vals))) {
      grid <- 0
    } else {
      rng <- range(vals, na.rm = TRUE)
      if (!all(is.finite(rng))) {
        grid <- 0
      } else if (identical(rng[1], rng[2])) {
        grid <- rep(rng[1], 2L) # duplicate to draw a flat line
      } else {
        grid <- seq(rng[1], rng[2], length.out = as.integer(grid_n))
      }
    }
    nd <- base[rep(1, length(grid)), , drop = FALSE]
    nd[[var]] <- grid
    nd$x <- nd[[var]]
    pred <- predict_wald(mod, nd)
    nd$fit <- pred$fit; nd$lwr <- pred$lwr; nd$upr <- pred$upr
    ggplot(nd, aes(x = x, y = fit)) +
      geom_line() +
      geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
      theme_minimal() +
      labs(title = as.character(mod_info$species), x = paste0(var, " (z)"), y = "Predicted events (offset = mean trap-nights)")
  }
}

predict_wald <- function(mod, newdata) {
  # Linear predictor and SE via model matrix and vcov; include offset from newdata
  linkinv <- mod$family$linkinv
  Xeta <- mm_and_eta(mod, newdata)
  V <- stats::vcov(mod)
  se_eta <- sqrt(pmax(0, rowSums((Xeta$mm %*% V) * Xeta$mm)))
  z <- stats::qnorm(0.975)
  fit <- linkinv(Xeta$eta)
  lwr <- linkinv(Xeta$eta - z * se_eta)
  upr <- linkinv(Xeta$eta + z * se_eta)
  list(fit = as.numeric(fit), lwr = as.numeric(lwr), upr = as.numeric(upr))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
