#' Ensemble pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param methods Methods to combine.
#' @param combine Combination strategy: mean or median.
#' @param auc_max_rank Fraction of top ranked genes for auc method.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_ensemble_matrix <- function(
  expr,
  genesets,
  methods = c("rank", "auc", "zmean"),
  combine = c("mean", "median"),
  auc_max_rank = 0.05,
  verbose = TRUE
) {
  methods <- unique(methods)
  combine <- match.arg(combine)

  mats <- list()
  for (m in methods) {
    mats[[m]] <- switch(
      m,
      rank = score_rank_matrix(expr, genesets, verbose = verbose),
      auc = score_auc_matrix(expr, genesets, auc_max_rank = auc_max_rank, verbose = verbose),
      zmean = score_zmean_matrix(expr, genesets, verbose = verbose),
      stop("Unsupported ensemble method: ", m, call. = FALSE)
    )
  }

  arr <- simplify2array(mats)
  if (length(dim(arr)) != 3L) {
    return(mats[[1]])
  }

  out <- if (combine == "mean") {
    apply(arr, c(1, 2), function(v) mean(v, na.rm = TRUE))
  } else {
    apply(arr, c(1, 2), function(v) stats::median(v, na.rm = TRUE))
  }

  rownames(out) <- rownames(mats[[1]])
  colnames(out) <- colnames(mats[[1]])
  out
}
