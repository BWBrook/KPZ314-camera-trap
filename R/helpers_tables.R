# R/helpers_tables.R
import::here(kable, .from = "knitr")
import::here(kable_styling, .from = "kableExtra")

# Simple kable wrapper with consistent styling
make_kable <- function(df, caption = NULL) {
  kable(df, caption = caption, digits = 2) |>
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
}
