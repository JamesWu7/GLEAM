#' Ensemble pathway scoring on matrix input
#'
#' @param expr Gene-by-cell matrix.
#' @param genesets Named list of pathways.
#' @param methods Methods to combine.
#' @param combine Combination strategy: mean or median.
#' @param standardize Ensemble harmonization (`zscore` or `rank`).
#' @param weights Optional named numeric vector of per-method weights.
#' @param verbose Whether to print progress.
#'
#' @return Pathway-by-cell numeric matrix.
#' @keywords internal
score_ensemble_matrix <- function(
  expr,
  genesets,
  methods = c("rank", "zscore", "mean"),
  combine = c("mean", "median"),
  standardize = c("zscore", "rank"),
  weights = NULL,
  verbose = TRUE
) {
  methods <- canonicalize_scoring_methods(methods)
  combine <- match.arg(combine)
  standardize <- match.arg(standardize)

  mats <- list()
  for (m in methods) {
    mats[[m]] <- switch(
      m,
      rank = score_rank_matrix(expr, genesets, verbose = verbose),
      mean = score_mean_matrix(expr, genesets, verbose = verbose),
      zscore = score_zscore_matrix(expr, genesets, verbose = verbose),
      scaled_mean = score_scaled_mean_matrix(expr, genesets, verbose = verbose),
      robust_mean = score_robust_matrix(expr, genesets, verbose = verbose),
      stop("Unsupported ensemble method: ", m, call. = FALSE)
    )
  }

  if (length(mats) == 1L) {
    return(mats[[1]])
  }

  mats <- lapply(mats, function(x) .harmonize_score_matrix(x, mode = standardize))
  w <- .resolve_ensemble_weights(methods = names(mats), weights = weights)

  arr <- simplify2array(mats)
  out <- if (combine == "mean" || !is.null(w)) {
    apply(arr, c(1, 2), function(v) {
      ok <- !is.na(v)
      if (!any(ok)) return(NA_real_)
      ww <- w[ok]
      if (all(ww == 0)) return(NA_real_)
      sum(v[ok] * ww) / sum(ww)
    })
  } else {
    apply(arr, c(1, 2), function(v) stats::median(v, na.rm = TRUE))
  }

  rownames(out) <- rownames(mats[[1]])
  colnames(out) <- colnames(mats[[1]])
  out
}

#' @keywords internal
.harmonize_score_matrix <- function(mat, mode = c("zscore", "rank")) {
  mode <- match.arg(mode)
  mat <- as.matrix(mat)
  out <- mat

  if (mode == "zscore") {
    for (i in seq_len(nrow(mat))) {
      v <- mat[i, ]
      mu <- mean(v, na.rm = TRUE)
      sdv <- stats::sd(v, na.rm = TRUE)
      if (is.na(sdv) || sdv == 0) {
        out[i, ] <- 0
      } else {
        out[i, ] <- (v - mu) / sdv
      }
    }
    return(out)
  }

  for (i in seq_len(nrow(mat))) {
    v <- mat[i, ]
    ok <- !is.na(v)
    if (!any(ok)) {
      out[i, ] <- NA_real_
      next
    }
    rv <- rep(NA_real_, length(v))
    rv[ok] <- rank(v[ok], ties.method = "average")
    n <- sum(ok)
    if (n > 1L) {
      rv[ok] <- (rv[ok] - 1) / (n - 1)
    } else {
      rv[ok] <- 0.5
    }
    out[i, ] <- rv
  }
  out
}

#' @keywords internal
.resolve_ensemble_weights <- function(methods, weights = NULL) {
  if (is.null(weights)) {
    out <- rep(1, length(methods))
    names(out) <- methods
    return(out)
  }
  if (!is.numeric(weights) || is.null(names(weights))) {
    stop("`ensemble_weights` must be a named numeric vector.", call. = FALSE)
  }
  if (any(weights < 0, na.rm = TRUE)) {
    stop("`ensemble_weights` cannot contain negative values.", call. = FALSE)
  }
  out <- rep(1, length(methods))
  names(out) <- methods
  shared <- intersect(methods, names(weights))
  out[shared] <- as.numeric(weights[shared])
  out
}
