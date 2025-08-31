# R/helpers_io.R
# read_mewc()  ── expert-checked MEWC table -------------------------------
# Expects a CSV with at least:
#   site, common, count, datetime, event, …
import::here(read_csv, col_factor, col_datetime, col_integer, 
             .from = "readr")
import::here(here, .from = "here")

#' Read expert-checked MEWC table
#' @param path Relative path to csv (default = data/mewc_master.csv)
#' @return tibble
#' @export
read_mewc <- function(path = here("data", "kpz314_2025_cam_data.csv")) {
  read_csv(
    path,
    show_col_types = FALSE,
    col_types = list(
      site = col_factor(),
      common  = col_factor(),
      count       = col_integer(),
      datetime   = col_datetime("%d/%m/%Y %H:%M:%S"),
      event       = col_integer()
    )
  )
}

# read_taxa()  ── species trait lookup ------------------------------------
# Expects at least: common, body_mass, trophic_group, …
#' @param path Relative path to csv (default = data/taxa.csv)
#' @export
read_taxa <- function(path = here("data", "kpz314_2025_species_list.csv")) {
  read_csv(path, show_col_types = FALSE)
}
