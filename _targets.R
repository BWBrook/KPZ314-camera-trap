# _targets.R  — executed only by staff *before* the lab -----------------
import::from("targets",   tar_option_set, tar_target, tar_source)
import::from("dplyr",     left_join, group_by, summarise, n, n_distinct, arrange, bind_rows)
import::from("readr",     read_csv)
import::from("here",      here)
import::from("tibble",    column_to_rownames)
import::from("tidyr",     pivot_wider, pivot_longer)
import::from("vegan",     vegdist)
import::from("ggplot2",   ggplot, aes, geom_col, scale_colour_manual)

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
  tar_target(gamma_by_habitat, gamma_hill(joined, group = "type", n_boot = 1000L, seed = 1L)),
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
  tar_target(nmds_stats,      bind_rows(nmds_bc$stats, nmds_jac$stats)),
  tar_target(beta_tests,      bind_rows(
                                    permanova_and_dispersion(beta_bc, site, group_col = "type", n_perm = 999L, seed = 1L),
                                    permanova_and_dispersion(beta_jac, site, group_col = "type", n_perm = 999L, seed = 1L)
                                  )),
  tar_target(beta_permanova,  beta_tests),
  tar_target(sp_table,        species_trend(joined,
                                  c("tasmanian_pademelon", "brushtail_possum"))),
  tar_target(fig_sp_trends, {
        dat <- pivot_longer(sp_table, -site)
        dat <- left_join(dat, site[, c("site", "type")], by = "site")
        ggplot(dat, aes(site, value, fill = name)) +
          geom_col(position = "dodge") +
          ggplot2::facet_wrap(~ type, ncol = 1)
      })
)
