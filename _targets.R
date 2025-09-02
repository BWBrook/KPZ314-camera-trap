# _targets.R  — executed only by staff *before* the lab -----------------
import::from("targets",   tar_option_set, tar_target, tar_source)
import::from("dplyr",     left_join, group_by, summarise, n, n_distinct, arrange, bind_rows)
import::from("readr",     read_csv)
import::from("here",      here)
import::from("tibble",    column_to_rownames)
import::from("tidyr",     pivot_wider, pivot_longer)
import::from("vegan",     vegdist)
import::from("ggplot2",   ggplot, aes, geom_col, scale_colour_manual)
import::from("tibble",    tibble)

tar_source("R") # load custom functions
tar_option_set(seed = 1L)

list(
  tar_target(site,            read_csv(here("data", "kpz314_2025_site_data.csv"),
                                  show_col_types = FALSE)),
  tar_target(raw,             read_mewc()),
  tar_target(taxa,            read_taxa()),
  # events + effort + RAI
  tar_target(events_default,  build_events(raw, gap_min = 5L)),
  tar_target(effort_site,     summarise_effort(site)),
  tar_target(rai,             compute_rai(events_default, effort_site)),
  tar_target(gap_values,      c(1L, 5L, 10L, 30L)),
  tar_target(events_gap,      build_events(raw, gap_min = gap_values),
                                  pattern = map(gap_values), iteration = "list"),
  tar_target(gap_sensitivity, summarise_gap_sensitivity(events_gap, effort_site)),
  tar_target(fig_gap,         gap_sensitivity$fig),
  tar_target(tab_effort,      make_kable(effort_site, caption = "Site effort manifest")),

  # join with site metadata for downstream analyses
  tar_target(joined,          left_join(events_default, site, by = "site")),
  tar_target(alpha,           calc_alpha(events_default)),
  tar_target(alpha_ci,        bootstrap_alpha(events_default, n_boot = 500L, seed = 1L)),
  tar_target(gamma_all,       n_distinct(events_default$common)),
  tar_target(gamma_region,    joined |> group_by(type) |>
                                  summarise(gamma = n_distinct(common))),
  tar_target(comm,            build_comm_matrix(events_default)),
  tar_target(fig_site_map,    plot_site_map(site)),
  tar_target(
    comm_region,
    joined |>
      group_by(type, common, event) |>
      summarise(n = 1L, .groups = "drop") |>
      group_by(type, common) |>
      summarise(events = n(), .groups = "drop") |>
      pivot_wider(
        names_from  = common,
        values_from = events,
        values_fill = 0L
      ) |>
      column_to_rownames("type") |>
      as.matrix()
  ),
  tar_target(rare_raw,        rarefaction_site(comm, step = 2L)),
  tar_target(rare_region,     rarefaction_by_group(comm, site, group_col = "type", step = 5L)),
  tar_target(gamma_by_habitat,
    {
      events_with_type <- dplyr::left_join(events_default, dplyr::select(site, site, type), by = "site")
      gamma_hill(events_with_type, group = "type", n_boot = 1000L, seed = 1L)
    }
  ),
  tar_target(gamma_all_ci,     gamma_hill(events_default, group = "all", n_boot = 1000L, seed = 1L)),
  # Distances and turnover visuals (paired: Bray–Curtis, Jaccard)
  tar_target(beta_bc,         compute_dist(comm, method = "bray")),
  tar_target(beta_jac,        compute_dist(comm, method = "jaccard")),
  tar_target(fig_heat_bc,     plot_heat_by_group(beta_bc, site, group_col = "type")),
  tar_target(fig_heat_jac,    plot_heat_by_group(beta_jac, site, group_col = "type")),
  # keep legacy name mapping to Bray–Curtis heatmap for backward compatibility
  tar_target(fig_heat,        fig_heat_bc),
  tar_target(nmds_bc,         nmds_with_hulls(beta_bc, site, group_col = "type", k = 2L, seed = 1L)),
  tar_target(nmds_jac,        nmds_with_hulls(beta_jac, site, group_col = "type", k = 2L, seed = 1L)),
  tar_target(fig_nmds,        nmds_bc$fig),           # preserve existing name for Bray–Curtis
  tar_target(fig_nmds_jac,    nmds_jac$fig),
  tar_target(fig_nmds_ellipse_bc,  nmds_with_ellipses(beta_bc, site, group_col = "type", k = 2L, seed = 1L)$fig),
  tar_target(nmds_stats,      bind_rows(nmds_bc$stats, nmds_jac$stats)),
  tar_target(beta_tests,      bind_rows(
                                    permanova_and_dispersion(beta_bc, site, group_col = "type", n_perm = 999L, seed = 1L),
                                    permanova_and_dispersion(beta_jac, site, group_col = "type", n_perm = 999L, seed = 1L)
                                  )),
  tar_target(beta_permanova,  beta_tests),
  tar_target(beta_r2_delta, {
      bt <- beta_tests
      bray <- bt[bt$distance == "bray", , drop = FALSE]
      jacc <- bt[bt$distance == "jaccard", , drop = FALSE]
      tibble::tibble(delta_r2 = bray$R2 - jacc$R2)
  }),
  tar_target(sp_table,        species_trend(joined,
                                  c("tasmanian_pademelon", "brushtail_possum"))),
  tar_target(fig_sp_trends, {
        dat <- pivot_longer(sp_table, -site)
        dat <- left_join(dat, site[, c("site", "type")], by = "site")
        ggplot(dat, aes(site, value, fill = name)) +
          geom_col(position = "dodge") +
          ggplot2::facet_wrap(~ type, ncol = 1)
      })

  ,
  # Habitat selection models (effort-offset GLMs)
  # Re-enable 'cat'; remove 'cwd' as predictor (only type, leaf_litter, can_foliage)
  tar_target(focal_species, c(
      "bennetts_wallaby","brushtail_possum","long_nosed_potoroo",
      "brush_bronzewing","cat"
    )),
  tar_target(glm_covars, c("type","leaf_litter","can_foliage")),
  tar_target(glm_data, prepare_species_glm_data(
      events_default, site, effort_site, species = focal_species, covars = glm_covars
    ), pattern = map(focal_species), iteration = "list"),
  tar_target(glm_models, fit_species_glm(
      df = glm_data, species = attr(glm_data, "species")
    ), pattern = map(glm_data), iteration = "list"),
  tar_target(glm_irrs, tidy_irrs(glm_models), pattern = map(glm_models), iteration = "list"),
  tar_target(glm_diag, diagnostics_table(glm_models), pattern = map(glm_models), iteration = "list"),
  tar_target(glm_irrs_all, bind_rows(glm_irrs)),
  tar_target(glm_diag_all, bind_rows(glm_diag)),
  tar_target(fig_glm_coef, plot_glm_coefs(glm_irrs), pattern = map(glm_irrs), iteration = "list"),
  tar_target(fig_glm_pd_type, partial_dependence(glm_models, df = NULL, var = "type"),
             pattern = map(glm_models), iteration = "list"),
  tar_target(fig_glm_pd_litter, partial_dependence(glm_models, df = NULL, var = "leaf_litter_z"),
             pattern = map(glm_models), iteration = "list"),
  tar_target(fig_glm_pd_cwd, partial_dependence(glm_models, df = NULL, var = "cwd_z"),
             pattern = map(glm_models), iteration = "list"),
  tar_target(fig_glm_pd_canopy, partial_dependence(glm_models, df = NULL, var = "can_foliage_z"),
             pattern = map(glm_models), iteration = "list")
  ,
  # Detectability: daily detection histories, ψ̂(t) curves, occupancy demo, and unknown sensitivity
  tar_target(detect_species, focal_species),
  tar_target(det_days_max, {
      md <- suppressWarnings(floor(min(site$op_days, na.rm = TRUE)))
      as.integer(max(1L, min(21L, ifelse(is.finite(md), md, 21L))))
    }
  ),
  tar_target(det_histories,
    detection_history_for_species(events_default, site, species = detect_species, max_days = det_days_max),
    pattern = map(detect_species), iteration = "list"
  ),
  tar_target(fig_det_heat,
    plot_detection_heatmap(det_histories$mat, det_histories$meta, det_histories$species),
    pattern = map(det_histories), iteration = "list"
  ),
  tar_target(psi_curves,
    naive_psi_curves(det_histories$mat, det_histories$meta),
    pattern = map(det_histories), iteration = "list"
  ),
  tar_target(psi_curves_all, bind_rows(psi_curves)),
  tar_target(fig_psi_curves,
    plot_psi_curves(psi_curves),
    pattern = map(psi_curves), iteration = "list"
  ),
  # Unknown-animal sensitivity
  tar_target(events_no_unknown, events_without_unknown(events_default)),
  tar_target(alpha_no_unknown,  calc_alpha(events_no_unknown)),
  tar_target(alpha_delta_unknown_tbl, alpha_delta_unknown(alpha, alpha_no_unknown)),
  tar_target(comm_no_unknown,    build_comm_matrix(events_no_unknown)),
  tar_target(beta_bc_no_unknown, compute_dist(comm_no_unknown, method = "bray")),
  tar_target(beta_jac_no_unknown, compute_dist(comm_no_unknown, method = "jaccard")),
  tar_target(beta_delta_unknown,  beta_delta_summary(comm, comm_no_unknown)),
  # Detection-rate ranking to auto-pick a pedagogical species
  tar_target(det_detect_summary,
    detection_rate_summary(det_histories$mat, det_histories$meta),
    pattern = map(det_histories), iteration = "list"
  ),
  tar_target(det_detect_summary_all, bind_rows(det_detect_summary)),
  tar_target(occ_species_auto, choose_occ_species(det_detect_summary_all)),

  # Occupancy demo: auto-picked species, ψ ~ type, p ~ 1
  tar_target(det_hist_occ,
    detection_history_for_species(events_default, site, species = occ_species_auto, max_days = det_days_max)
  ),
  tar_target(umf_occ,      occu_build_umf(det_hist_occ$mat, det_hist_occ$meta)),
  tar_target(occu_demo,    fit_simple_occu(umf_occ, with_type = TRUE)),
  tar_target(tab_occu,     occu_demo$tidy),
  tar_target(fig_occu,     occu_demo$fig),
  # Inspection targets for occupancy stability
  tar_target(occ_hist_long, detection_history_long(det_hist_occ$mat, det_hist_occ$meta)),
  tar_target(occ_site_summary, summarise_occ_site(det_hist_occ$mat, det_hist_occ$meta))
)
