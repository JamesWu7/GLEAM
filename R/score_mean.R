#' Mean pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_mean_matrix <- function(expr, genesets, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  if (verbose) message("[GLEAM] scoring mean method...")

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(expr), nomatch = 0L)
    idx <- idx[idx > 0L]
    col_means_by_idx(expr, idx)
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(expr)
  mat
}
