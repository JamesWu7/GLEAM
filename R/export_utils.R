#' Collect score objects from multiple methods
#'
#' @param ... Named `gleam_score` objects.
#'
#' @return Named list of score objects.
#' @export
collect_scores <- function(...) {
  xs <- list(...)
  if (length(xs) < 1L) stop("Provide at least one score object.", call. = FALSE)
  ok <- vapply(xs, function(x) inherits(x, "gleam_score") || inherits(x, "scpathway_score"), logical(1))
  if (!all(ok)) stop("All inputs must be score objects returned by score_pathway().", call. = FALSE)
  xs
}

#' Pivot score matrix to long format
#'
#' @param score `gleam_score` object.
#'
#' @return Long-format data.frame with `cell_id`, `pathway`, and `score`.
#' @export
pivot_scores_long <- function(score) {
  check_score_object(score)
  mat <- score$score
  expand <- expand.grid(pathway = rownames(mat), cell_id = colnames(mat), stringsAsFactors = FALSE)
  expand$score <- mapply(function(pw, cid) mat[pw, cid], expand$pathway, expand$cell_id)
  expand
}

#' Summarize score table
#'
#' @param score `gleam_score` object.
#' @param by Grouping columns.
#' @param fun Aggregation function.
#' @param threshold Threshold used for fraction.
#'
#' @return Summarized data.frame.
#' @export
summarize_scores <- function(score, by = c("sample", "group"), fun = c("mean", "median", "fraction"), threshold = 0) {
  fun <- match.arg(fun)
  aggregate_pathway(score = score, by = by, fun = fun, threshold = threshold, long = TRUE)
}

#' Export score table
#'
#' @param x `gleam_score` object or data.frame.
#' @param file Output file path.
#' @param format Output format: `csv` or `tsv`.
#' @param include_meta Whether to include metadata columns when input is score object.
#'
#' @return Invisibly the exported data.frame.
#' @export
export_scores <- function(x, file, format = c("csv", "tsv"), include_meta = TRUE) {
  format <- match.arg(format)
  df <- if (inherits(x, "gleam_score") || inherits(x, "scpathway_score")) {
    long <- pivot_scores_long(x)
    if (include_meta) {
      meta <- x$meta
      meta$cell_id <- rownames(meta)
      merge(long, meta, by = "cell_id", all.x = TRUE)
    } else {
      long
    }
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("`x` must be a score object or data.frame.", call. = FALSE)
  }

  if (format == "csv") {
    utils::write.csv(df, file = file, row.names = FALSE)
  } else {
    utils::write.table(df, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  invisible(df)
}

#' Compare multiple scoring methods
#'
#' @param ... Named score objects from different methods.
#' @param pathway Optional pathway subset.
#' @param summary_fun Summary function (`mean` or `median`).
#'
#' @return Data.frame of per-method pathway summaries.
#' @export
compare_scoring_methods <- function(..., pathway = NULL, summary_fun = c("mean", "median")) {
  xs <- collect_scores(...)
  summary_fun <- match.arg(summary_fun)

  rows <- lapply(names(xs), function(nm) {
    s <- xs[[nm]]
    mat <- s$score
    if (!is.null(pathway)) {
      keep <- intersect(pathway, rownames(mat))
      mat <- mat[keep, , drop = FALSE]
    }
    val <- if (summary_fun == "mean") {
      rowMeans(mat, na.rm = TRUE)
    } else {
      apply(mat, 1, stats::median, na.rm = TRUE)
    }
    data.frame(method = nm %||% s$method, pathway = names(val), value = as.numeric(val), stringsAsFactors = FALSE)
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
