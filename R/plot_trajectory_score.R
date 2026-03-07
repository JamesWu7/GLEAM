#' Plot pathway score on trajectory embedding
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param embeddings Embedding matrix with at least 2 columns.
#' @param reduction Reduction name for Seurat extraction when `embeddings` is NULL.
#' @param object Optional Seurat object for reduction extraction.
#'
#' @return A `ggplot` object.
#' @export
plot_trajectory_score <- function(score, pathway, embeddings = NULL, reduction = "umap", object = NULL) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  emb <- if (!is.null(embeddings)) {
    as.matrix(embeddings)
  } else {
    if (is.null(object)) stop("Provide `embeddings` or a Seurat `object`.", call. = FALSE)
    extract_reduction(object = object, reduction = reduction)
  }

  if (ncol(emb) < 2L) stop("`embeddings` must have at least 2 columns.", call. = FALSE)
  emb <- emb[colnames(score$score), , drop = FALSE]

  df <- data.frame(
    dim1 = emb[, 1],
    dim2 = emb[, 2],
    score = as.numeric(score$score[pathway, ]),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, color = score)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#3b4cc0", mid = "#f7f7f7", high = "#b40426", midpoint = 0) +
    ggplot2::labs(title = paste("Trajectory map:", pathway), x = colnames(emb)[1], y = colnames(emb)[2]) +
    .theme_gleam()
}
