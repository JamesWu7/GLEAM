#' @keywords internal
as_expr_matrix <- function(expr) {
  if (inherits(expr, "dgCMatrix")) {
    return(expr)
  }
  if (is.matrix(expr)) {
    return(expr)
  }
  stop("`expr` must be a matrix or Matrix::dgCMatrix.", call. = FALSE)
}

#' @keywords internal
as_dense_matrix <- function(x) {
  if (inherits(x, "dgCMatrix")) {
    return(as.matrix(x))
  }
  x
}

#' @keywords internal
col_means_by_idx <- function(expr, idx) {
  if (length(idx) == 0L) {
    return(rep(NA_real_, ncol(expr)))
  }
  if (inherits(expr, "dgCMatrix")) {
    return(as.numeric(Matrix::colMeans(expr[idx, , drop = FALSE], na.rm = TRUE)))
  }
  as.numeric(colMeans(expr[idx, , drop = FALSE], na.rm = TRUE))
}

#' @keywords internal
col_ranks_desc <- function(expr) {
  dense <- as_dense_matrix(expr)
  out <- apply(dense, 2, function(v) rank(-v, ties.method = "average"))
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1L)
  }
  rownames(out) <- rownames(dense)
  colnames(out) <- colnames(dense)
  out
}

#' @keywords internal
zscore_by_row <- function(expr) {
  dense <- as_dense_matrix(expr)
  mu <- rowMeans(dense, na.rm = TRUE)
  sdv <- apply(dense, 1, stats::sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  (dense - mu) / sdv
}
