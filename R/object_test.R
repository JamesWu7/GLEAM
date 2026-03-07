#' Construct `scpathway_test` object
#'
#' @param table Differential pathway result table.
#' @param level Comparison level.
#' @param method Statistical method.
#' @param comparison Comparison details.
#' @param params Parameter list.
#'
#' @return Object of class `scpathway_test`.
#' @keywords internal
new_scpathway_test <- function(table, level, method, comparison = list(), params = list()) {
  structure(
    list(
      table = table,
      level = level,
      method = method,
      comparison = comparison,
      params = params
    ),
    class = "scpathway_test"
  )
}
