# R/helpers_text.R — small text helpers for practicum captions

#' Turnover caption explaining distance choices and dispersion caveat
#'
#' @param distance Optional character: "bray" or "jaccard" to tailor wording
#' @return Length-1 character vector suitable for ggplot captions or inline text
#' @export
turnover_caption <- function(distance = NULL) {
  base <- "Bray–Curtis (abundance-sensitive) highlights dominance shifts; Jaccard (presence–absence) emphasises composition regardless of counts. Report PERMANOVA R² and permutation p for each distance. If dispersion p < 0.05, interpret the location effect cautiously (groups differ in spread)."
  if (is.null(distance)) return(base)
  distance <- tolower(as.character(distance))
  if (distance == "bray") {
    return("Bray–Curtis reflects abundance and dominance differences between sites.")
  } else if (distance == "jaccard") {
    return("Jaccard reflects presence–absence (composition) irrespective of abundances.")
  }
  base
}

