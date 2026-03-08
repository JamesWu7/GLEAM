#' Aggregate pathway scores
#'
#' @param score `gleam_score` object.
#' @param by Grouping variable(s): metadata column name(s), or vector.
#' @param fun Aggregation function: mean, median, fraction, sum, or a function.
#' @param threshold Threshold used by `fraction`.
#' @param long Return long-format table if TRUE.
#'
#' @return Aggregated table (long) or pathway-by-group matrix (wide).
#' @keywords internal
aggregate_pathway <- function(score, by, fun = c("mean", "median", "fraction", "sum"), threshold = 0, long = TRUE) {
  check_score_object(score)
  fun_name <- NULL
  if (is.function(fun)) {
    fun_name <- "custom"
  } else {
    fun <- match.arg(fun)
    fun_name <- fun
  }

  meta <- score$meta
  mat <- score$score

  if (is.character(by)) {
    check_required_columns(meta, by)
    grp_df <- meta[, by, drop = FALSE]
  } else if (is.vector(by) && length(by) == ncol(mat)) {
    grp_df <- data.frame(group = by, stringsAsFactors = FALSE, row.names = rownames(meta))
  } else {
    stop("`by` must be metadata column name(s) or a vector of length ncol(score).", call. = FALSE)
  }

  grp_key <- interaction(grp_df, drop = TRUE, lex.order = TRUE)
  grp_levels <- levels(grp_key)

  agg_fun <- if (is.function(fun)) {
    function(x) apply(x, 1, fun)
  } else {
    switch(
      fun,
      mean = function(x) rowMeans(x, na.rm = TRUE),
      median = function(x) apply(x, 1, stats::median, na.rm = TRUE),
      fraction = function(x) rowMeans(x > threshold, na.rm = TRUE),
      sum = function(x) rowSums(x, na.rm = TRUE)
    )
  }

  wide <- sapply(grp_levels, function(gl) {
    idx <- grp_key == gl
    agg_fun(mat[, idx, drop = FALSE])
  })
  if (is.null(dim(wide))) {
    wide <- matrix(wide, ncol = 1L)
  }
  rownames(wide) <- rownames(mat)
  colnames(wide) <- grp_levels

  if (!long) {
    return(wide)
  }

  group_map <- unique(data.frame(group_key = as.character(grp_key), grp_df, stringsAsFactors = FALSE))

  out <- do.call(rbind, lapply(seq_len(ncol(wide)), function(i) {
    key <- colnames(wide)[i]
    ginfo <- group_map[group_map$group_key == key, , drop = FALSE]
    data.frame(
      pathway = rownames(wide),
      group_key = key,
      value = as.numeric(wide[, i]),
      ginfo[rep(1L, nrow(wide)), setdiff(colnames(ginfo), "group_key"), drop = FALSE],
      stringsAsFactors = FALSE
    )
  }))

  rownames(out) <- NULL
  out$summary_method <- fun_name
  out
}
