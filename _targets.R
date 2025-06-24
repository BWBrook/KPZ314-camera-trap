# _targets.R  — executed only by staff *before* the lab -----------------
import::from("targets", tar_option_set, tar_target)
import::here(read_mewc, read_taxa, .from = "R/helpers_io.R")
import::here(calc_alpha, build_comm_matrix, .from = "R/helpers_metrics.R")
import::here(plot_turnover_heatmap, plot_nmds, .from = "R/helpers_plots.R")
import::here(calc_rarefaction, permanova_region, species_trend, .from = "R/helpers_extend.R")

tar_option_set(packages = c("dplyr", "tidyr", "vegan", "ggplot2",
                            "here", "reshape2"))

list(
  tar_target(raw,      read_mewc()),
  tar_target(taxa,     read_taxa()),
  tar_target(site,     readr::read_csv(here::here("data",
                                                  "site_descriptor.csv"),
                                       show_col_types = FALSE)),
  tar_target(joined,   dplyr::left_join(raw, site, by = "camera_site")),
  tar_target(alpha,    calc_alpha(joined)),
  tar_target(gamma_all, dplyr::n_distinct(raw$class_name)),
  tar_target(gamma_region, joined |> dplyr::group_by(region) |>
             dplyr::summarise(gamma = dplyr::n_distinct(class_name))),
  tar_target(comm,     build_comm_matrix(joined)),
  tar_target(
    comm_region,
    joined |>
      dplyr::group_by(region, class_name, event) |>
      dplyr::summarise(n = 1L, .groups = "drop") |>
      dplyr::group_by(region, class_name) |>
      dplyr::summarise(events = n(), .groups = "drop") |>
      tidyr::pivot_wider(
        names_from  = class_name,
        values_from = events,
        values_fill = 0L
      ) |>
      tibble::column_to_rownames("region") |>
      as.matrix()
  ),
  tar_target(rare_raw,    calc_rarefaction(comm)),
  tar_target(rare_region, calc_rarefaction(comm_region)),
  tar_target(beta_permanova, permanova_region(beta_bc, site)),           
  tar_target(beta_bc,  vegan::vegdist(comm, "bray")),
  tar_target(fig_heat, plot_turnover_heatmap(beta_bc, site)),
  tar_target(fig_nmds, plot_nmds(beta_bc, joined)),
  tar_target(sp_table, species_trend(joined,
                     c("tasmanian_pademelon", "brushtail_possum"))),
  tar_target(fig_sp_trends, {
  ggplot2::ggplot(tidyr::pivot_longer(sp_table, -camera_site),
           ggplot2::aes(camera_site, value, fill = name)) + 
           ggplot2::geom_col(position = "dodge")})
)
