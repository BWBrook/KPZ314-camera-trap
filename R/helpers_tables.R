# R/helpers_tables.R
# Simple kable wrapper with consistent styling
make_kable <- function(df, caption = NULL) {
  import::from("knitr", kable)
  import::from("kableExtra", kable_styling)

  kable(df, caption = caption, digits = 2) |>
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
}
