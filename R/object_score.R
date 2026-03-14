#' Construct score object
#'
#' @param score Pathway-by-cell matrix.
#' @param meta Cell metadata `data.frame` aligned to score columns.
#' @param method Scoring method.
#' @param geneset_name Geneset label.
#' @param geneset_info Geneset metadata.
#' @param params Parameter list.
#'
#' @return Object of class `gleam_score` and `gleam_score`.
#' @keywords internal
new_gleam_score <- function(score, meta, method, geneset_name, geneset_info, params = list()) {
  if (is.null(colnames(score))) {
    stop("`score` must have cell names in colnames().", call. = FALSE)
  }
  if (nrow(meta) != ncol(score)) {
    stop("`meta` row count must equal ncol(score).", call. = FALSE)
  }
  if (!identical(rownames(meta), colnames(score))) {
    meta <- meta[colnames(score), , drop = FALSE]
  }

  obj <- structure(
    list(
      score = score,
      meta = meta,
      method = method,
      geneset_name = geneset_name,
      geneset_info = geneset_info,
      params = params
    ),
    class = "gleam_score"
  )
  obj
}
