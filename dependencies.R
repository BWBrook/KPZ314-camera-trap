project_dependencies <- function() {
  c(
    # core tidy and plotting
    "dplyr","tidyr","tibble","ggplot2","reshape2",
    # analysis
    "vegan","sf","cli","rlang","withr","lubridate",
    # GLM NB
    "MASS",
    # occupancy demo
    "unmarked"
  )
}
