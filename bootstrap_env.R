if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
import::from(renv, init, install, snapshot, status)

required <- c(
  "targets", "import", "here", "tidyverse", "vegan", 
  "broom", "patchwork", "ggvegan", "plotly"
)

init(bare = TRUE) |> invisible()
missing <- setdiff(required, status()$library$Package)
if (length(missing)) install(missing)
snapshot()
