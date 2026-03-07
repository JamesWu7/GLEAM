#' Singscore-like pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_singscore_like_matrix <- function(expr, genesets, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  ranks <- col_ranks_desc(expr)
  ng <- nrow(ranks)
  if (verbose) message("[GLEAM] scoring singscore_like method...")

  out <- lapply(names(genesets), function(pw) {
    idx <- match(genesets[[pw]], rownames(ranks), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) return(rep(NA_real_, ncol(ranks)))
    raw <- colMeans(ranks[idx, , drop = FALSE], na.rm = TRUE)
    min_s <- (length(idx) + 1) / 2
    max_s <- ng - min_s
    scaled <- 1 - (raw - min_s) / pmax(1e-9, (max_s - min_s))
    as.numeric(scaled)
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(expr)
  mat
}

#' ssGSEA-like pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param alpha Rank-weight exponent.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_ssgsea_like_matrix <- function(expr, genesets, alpha = 0.25, verbose = TRUE) {
  expr <- as_expr_matrix(expr)
  dense <- as_dense_matrix(expr)
  if (verbose) message("[GLEAM] scoring ssgsea_like method...")

  ord_idx <- apply(dense, 2, order, decreasing = TRUE)
  if (is.null(dim(ord_idx))) ord_idx <- matrix(ord_idx, ncol = 1)

  out <- lapply(names(genesets), function(pw) {
    gs <- genesets[[pw]]
    v <- rep(NA_real_, ncol(dense))
    for (j in seq_len(ncol(dense))) {
      ord <- ord_idx[, j]
      genes_ord <- rownames(dense)[ord]
      hit <- genes_ord %in% gs
      if (!any(hit)) {
        v[j] <- NA_real_
        next
      }
      rank_pos <- seq_along(genes_ord)
      w_hit <- (length(genes_ord) - rank_pos + 1)^alpha
      phit <- cumsum(ifelse(hit, w_hit, 0)) / sum(w_hit[hit])
      pmiss <- cumsum(ifelse(!hit, 1, 0)) / sum(!hit)
      es <- phit - pmiss
      v[j] <- max(es, na.rm = TRUE)
    }
    v
  })

  mat <- do.call(rbind, out)
  rownames(mat) <- names(genesets)
  colnames(mat) <- colnames(dense)
  mat
}
