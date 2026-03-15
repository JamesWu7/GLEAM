#' Scaled mean pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_scaled_mean_matrix <- function(expr, genesets, verbose = TRUE) {
  if (verbose) message("[GLEAM] scoring scaled_mean method...")
  score_zscore_matrix(expr = expr, genesets = genesets, verbose = FALSE)
}

#' Robust normalized pathway scoring
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_robust_matrix <- function(expr, genesets, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  dense <- as_dense_matrix(expr)
  if (verbose) message("[GLEAM] scoring robust_mean method...")

  med <- apply(dense, 1, stats::median, na.rm = TRUE)
  madv <- apply(dense, 1, stats::mad, na.rm = TRUE)
  madv[madv == 0 | is.na(madv)] <- 1
  rz <- (dense - med) / madv

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(rz), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) return(rep(NA_real_, ncol(rz)))
    as.numeric(colMeans(rz[idx, , drop = FALSE], na.rm = TRUE))
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(dense)
  mat
}
