# _targets.R  — executed only by staff *before* the lab -----------------
import::from(targets, tar_option_set, tar_target)
import::here(read_mewc, read_taxa, .from = "R/helpers_io.R")
import::here(calc_alpha, build_comm_matrix, .from = "R/helpers_metrics.R")
import::here(plot_turnover_heatmap, plot_nmds, .from = "R/helpers_plots.R")

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
  tar_target(comm,     build_comm_matrix(joined)),
  tar_target(beta_bc,  vegan::vegdist(comm, "bray")),
  tar_target(fig_heat, plot_turnover_heatmap(beta_bc, site)),
  tar_target(fig_nmds, plot_nmds(beta_bc, joined))
)
