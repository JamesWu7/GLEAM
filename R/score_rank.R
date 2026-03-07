#' Rank-based pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param normalize Whether to normalize rank score to the 0-1 range.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_rank_matrix <- function(expr, genesets, normalize = TRUE, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  ranks <- col_ranks_desc(expr)
  ng <- nrow(expr)

  if (verbose) {
    message("[GLEAM] scoring rank method...")
  }

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(ranks), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) {
      return(rep(NA_real_, ncol(ranks)))
    }
    v <- colMeans(ranks[idx, , drop = FALSE], na.rm = TRUE)
    if (normalize && ng > 1L) {
      v <- 1 - (v - 1) / (ng - 1)
    }
    as.numeric(v)
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(expr)
  mat
}
