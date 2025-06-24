# R/helpers_io.R
# read_mewc()  ── expert-checked MEWC table -------------------------------
# Expects a CSV with at least:
#   camera_site, class_name, event, timestamp, …
import::from("readr", read_csv)
import::from("here", here)

#' Read expert-checked MEWC table
#' @param path Relative path to csv (default = data/mewc_master.csv)
#' @return tibble
#' @export
read_mewc <- function(path = here("data", "mewc_master.csv")) {
  read_csv(
    path,
    show_col_types = FALSE,
    col_types = list(
      camera_site = readr::col_factor(),
      class_name  = readr::col_factor(),
      event       = readr::col_integer(),
      timestamp   = readr::col_datetime("%d/%m/%Y %H:%M:%S")
    )
  )
}

# read_taxa()  ── species trait lookup ------------------------------------
# Expects at least: class_name, body_mass, trophic_group, …
#' @param path Relative path to csv (default = data/taxa.csv)
#' @export
read_taxa <- function(path = here("data", "taxa.csv")) {
  read_csv(path, show_col_types = FALSE)
}
