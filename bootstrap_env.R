if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
import::from("renv", init, install, snapshot, status)

required <- c(
  "targets", "import", "dplyr", "readr", "here", "tidyr", "vegan", "ggplot2", 
  "tibble", "reshape2", "quarto", "knitr", "kableExtra", "plotly", "ggspatial"
)

init(bare = TRUE) |> invisible()
missing <- setdiff(required, status()$library$Package)
if (length(missing)) install(missing)
snapshot()

quarto::check_newer_version()