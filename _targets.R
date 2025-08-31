# _targets.R  — executed only by staff *before* the lab -----------------
import::from("targets",   tar_option_set, tar_target, tar_source)
import::from("dplyr",     left_join, group_by, summarise, n, n_distinct, arrange)
import::from("readr",     read_csv)
import::from("here",      here)
import::from("tibble",    column_to_rownames)
import::from("tidyr",     pivot_wider, pivot_longer)
import::from("vegan",     vegdist)
import::from("ggplot2",   ggplot, aes, geom_col, scale_colour_manual)

tar_source("R") # load custom functions
tar_option_set(seed = 1L)

list(
  tar_target(site,     read_csv(here("data", "kpz314_2025_site_data.csv"),
                                show_col_types = FALSE)),
  tar_target(raw,      read_mewc()),
  tar_target(taxa,     read_taxa()),
  # events + effort + RAI
  tar_target(events_default, build_events(raw, gap_min = 5L)),
  tar_target(effort_site,    summarise_effort(site)),
  tar_target(rai,            compute_rai(events_default, effort_site)),
  tar_target(gap_values,     c(1L, 5L, 10L, 30L)),
  tar_target(events_gap,     build_events(raw, gap_min = gap_values),
             pattern = map(gap_values), iteration = "list"),
  tar_target(gap_sensitivity, summarise_gap_sensitivity(events_gap, effort_site)),
  tar_target(fig_gap,         gap_sensitivity$fig),
  tar_target(tab_effort,      make_kable(effort_site, caption = "Site effort manifest")),

  # join with site metadata for downstream analyses
  tar_target(joined,   left_join(events_default, site, by = "site")),
  tar_target(alpha,    calc_alpha(joined)),
  tar_target(gamma_all, n_distinct(events_default$common)),
  tar_target(gamma_region, joined |> group_by(type) |>
                summarise(gamma = n_distinct(common))),
  tar_target(comm,     build_comm_matrix(joined)),
  tar_target(fig_site_map, plot_site_map(site)),
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
  tar_target(rare_raw,    calc_rarefaction(comm)),
  tar_target(rare_region, {
    p <- calc_rarefaction(comm_region)
    p + ggplot2::aes(colour = site) +
      ggplot2::scale_colour_manual(
        values = c(dry = "sienna3", wet = "steelblue4"),
        name = "Type"
      )
  }),
  tar_target(beta_permanova, permanova_region(beta_bc, site)),           
  tar_target(beta_bc,  vegdist(comm, "bray")),
  tar_target(fig_heat, plot_turnover_heatmap(beta_bc, site)),
  tar_target(fig_nmds, plot_nmds(beta_bc, joined)),
  tar_target(sp_table, species_trend(joined,
             c("tasmanian_pademelon", "brushtail_possum"))),
  tar_target(fig_sp_trends, {
        dat <- pivot_longer(sp_table, -site)
        dat <- left_join(dat, site[, c("site", "type")], by = "site")
        ggplot(dat, aes(site, value, fill = name)) +
          geom_col(position = "dodge") +
          ggplot2::facet_wrap(~ type, ncol = 1)
      })
)
