if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", "import")
import::from("renv", init, install, snapshot, status)

required <- c(
  "targets", "import", "dplyr", "readr", "here", "tidyr", "vegan", "ggplot2", 
  "tibble", "reshape2", "knitr", "kableExtra"
)

init(bare = TRUE) |> invisible()
missing <- setdiff(required, status()$library$Package)
if (length(missing)) install(missing)
snapshot()
