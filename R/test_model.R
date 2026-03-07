#' Experimental linear model placeholder
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param formula Formula for model.
#' @param data Data frame.
#'
#' @return A placeholder list.
#' @keywords internal
test_pathway_lm_placeholder <- function(score_mat, formula, data) {
  warning("Linear model interface is experimental in v1 and not fully implemented.", call. = FALSE)
  list(score_mat = score_mat, formula = formula, data = data)
}
