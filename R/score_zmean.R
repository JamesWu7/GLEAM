#' Z-score pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_zmean_matrix <- function(expr, genesets, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  z <- zscore_by_row(expr)

  if (verbose) {
    message("[GLEAM] scoring zscore method...")
  }

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(z), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) {
      return(rep(NA_real_, ncol(z)))
    }
    as.numeric(colMeans(z[idx, , drop = FALSE], na.rm = TRUE))
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(expr)
  mat
}
