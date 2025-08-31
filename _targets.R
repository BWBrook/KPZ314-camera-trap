# _targets.R  — executed only by staff *before* the lab -----------------
import::from("targets", tar_option_set, tar_target)
import::from("dplyr", left_join, group_by, summarise, n, n_distinct)
import::from("readr", read_csv)
import::from("here", here)
import::from("tibble", column_to_rownames)
import::from("tidyr", pivot_wider, pivot_longer)
import::from("vegan", vegdist)
import::from("ggplot2", ggplot, aes, geom_col)

import::here(read_mewc, read_taxa, .from = "R/helpers_io.R")
import::here(calc_alpha, build_comm_matrix, .from = "R/helpers_metrics.R")
import::here(plot_turnover_heatmap, plot_nmds, .from = "R/helpers_plots.R")
import::here(calc_rarefaction, permanova_region, species_trend, 
             .from = "R/helpers_extend.R")

list(
  tar_target(raw,      read_mewc()),
  tar_target(taxa,     read_taxa()),
  tar_target(site,     read_csv(here("data", "kzp314_2025_site_data.csv"),
                                show_col_types = FALSE)),
  tar_target(joined,   left_join(raw, site, by = "site")),
  tar_target(alpha,    calc_alpha(joined)),
  tar_target(gamma_all, n_distinct(raw$common)),
  tar_target(gamma_region, joined |> group_by(region) |>
                summarise(gamma = n_distinct(common))),
  tar_target(comm,     build_comm_matrix(joined)),
  tar_target(
    comm_region,
    joined |>
      group_by(region, common, event) |>
      summarise(n = 1L, .groups = "drop") |>
      group_by(region, common) |>
      summarise(events = n(), .groups = "drop") |>
      pivot_wider(
        names_from  = common,
        values_from = events,
        values_fill = 0L
      ) |>
      column_to_rownames("region") |>
      as.matrix()
  ),
  tar_target(rare_raw,    calc_rarefaction(comm)),
  tar_target(rare_region, calc_rarefaction(comm_region)),
  tar_target(beta_permanova, permanova_region(beta_bc, site)),           
  tar_target(beta_bc,  vegdist(comm, "bray")),
  tar_target(fig_heat, plot_turnover_heatmap(beta_bc, site)),
  tar_target(fig_nmds, plot_nmds(beta_bc, joined)),
  tar_target(sp_table, species_trend(joined,
             c("tasmanian_pademelon", "brushtail_possum"))),
  tar_target(fig_sp_trends, {
        ggplot(pivot_longer(sp_table, -site),
        aes(site, value, fill = name)) + 
        geom_col(position = "dodge")})
)
