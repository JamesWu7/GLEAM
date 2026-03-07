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
  dense <- as_dense_matrix(expr)
  if (verbose) message("[GLEAM] scoring mean method...")

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(dense), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) return(rep(NA_real_, ncol(dense)))
    as.numeric(colMeans(dense[idx, , drop = FALSE], na.rm = TRUE))
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(dense)
  mat
}
