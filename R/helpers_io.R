# R/helpers_io.R
# read_mewc()  ── expert-checked MEWC table -------------------------------
# Expects a CSV with at least:
#   camera_site, class_name, event, timestamp, …
import::here(read_csv, col_factor, col_datetime, col_integer, 
             .from = "readr")
import::here(here, .from = "here")

#' Read expert-checked MEWC table
#' @param path Relative path to csv (default = data/mewc_master.csv)
#' @return tibble
#' @export
read_mewc <- function(path = here("data", "kpz314_2025_cam_dataset.csv")) {
  read_csv(
    path,
    show_col_types = FALSE,
    col_types = list(
      camera_site = col_factor(),
      class_name  = col_factor(),
      event       = col_integer(),
      timestamp   = col_datetime("%d/%m/%Y %H:%M:%S")
    )
  )
}

# read_taxa()  ── species trait lookup ------------------------------------
# Expects at least: class_name, body_mass, trophic_group, …
#' @param path Relative path to csv (default = data/taxa.csv)
#' @export
read_taxa <- function(path = here("data", "kzp314_2025_species_list.csv")) {
  read_csv(path, show_col_types = FALSE)
}
