#' @export
print.scpathway_score <- function(x, ...) {
  cat("<scpathway_score>\n")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  pathways: ", nrow(x$score), "\n", sep = "")
  cat("  cells: ", ncol(x$score), "\n", sep = "")
  invisible(x)
}

#' @export
print.scpathway_test <- function(x, ...) {
  cat("<scpathway_test>\n")
  cat("  level: ", x$level, "\n", sep = "")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  rows: ", nrow(x$table), "\n", sep = "")
  print(utils::head(x$table, 6))
  invisible(x)
}

#' @export
summary.scpathway_test <- function(object, ...) {
  tbl <- object$table
  out <- list(
    level = object$level,
    method = object$method,
    n_pathway_tests = nrow(tbl),
    n_significant = sum(tbl$p_adj < 0.05, na.rm = TRUE)
  )
  class(out) <- "summary.scpathway_test"
  out
}
