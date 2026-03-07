#' AUC-like pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param auc_max_rank Fraction of top-ranked genes considered active.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_auc_matrix <- function(expr, genesets, auc_max_rank = 0.05, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  ranks <- col_ranks_desc(expr)
  top_n <- max(1L, ceiling(nrow(expr) * auc_max_rank))

  if (verbose) {
    message(sprintf("[GLEAM] scoring auc method (top_n=%d)...", top_n))
  }

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(ranks), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) {
      return(rep(NA_real_, ncol(ranks)))
    }
    as.numeric(colMeans(ranks[idx, , drop = FALSE] <= top_n, na.rm = TRUE))
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(expr)
  mat
}
