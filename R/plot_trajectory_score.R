#' Plot pathway score on trajectory embedding
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param embeddings Embedding matrix with at least 2 columns.
#' @param reduction Reduction name for Seurat extraction when `embeddings` is NULL.
#' @param object Optional Seurat object for reduction extraction.
#' @param point_size Point size.
#' @param alpha Point alpha.
#' @param palette Continuous palette name or custom colors.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_trajectory_score <- function(
  score,
  pathway,
  embeddings = NULL,
  reduction = "umap",
  object = NULL,
  point_size = 1.1,
  alpha = 0.9,
  palette = "gleam_continuous",
  theme_params = list()
) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  emb <- if (!is.null(embeddings)) {
    as.matrix(embeddings)
  } else {
    if (is.null(object)) stop("Provide `embeddings` or a Seurat `object`.", call. = FALSE)
    extract_reduction(object = object, reduction = reduction)
  }

  if (ncol(emb) < 2L) stop("`embeddings` must have at least 2 columns.", call. = FALSE)
  if (is.null(rownames(emb))) {
    if (nrow(emb) != ncol(score$score)) {
      stop("Embedding has no rownames; nrow(embeddings) must equal ncol(score).", call. = FALSE)
    }
    rownames(emb) <- colnames(score$score)
  }
  if (!all(colnames(score$score) %in% rownames(emb))) {
    stop("Embedding rownames must include all score cell IDs.", call. = FALSE)
  }
  emb <- emb[colnames(score$score), , drop = FALSE]

  df <- data.frame(
    dim1 = emb[, 1],
    dim2 = emb[, 2],
    score = as.numeric(score$score[pathway, ]),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$dim1, y = ggplot2::.data$dim2, color = ggplot2::.data$score)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::labs(title = paste("Trajectory map:", pathway), x = colnames(emb)[1], y = colnames(emb)[2]) +
    do.call(gleam_theme, tp)
}
