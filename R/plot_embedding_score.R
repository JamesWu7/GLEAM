#' Plot pathway score on embedding coordinates
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param embedding Embedding matrix with at least 2 columns.
#' @param object Optional Seurat object.
#' @param reduction Reduction name used when `embedding` is NULL.
#' @param split.by Optional metadata column for faceting.
#' @param point_size Point size.
#' @param alpha Point alpha.
#' @param palette Continuous palette name or custom colors.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_embedding_score <- function(
  score,
  pathway,
  embedding = NULL,
  object = NULL,
  reduction = "umap",
  split.by = NULL,
  point_size = 1.1,
  alpha = 0.9,
  palette = "gleam_continuous",
  theme_params = list()
) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)
  emb <- if (!is.null(embedding)) {
    as.matrix(embedding)
  } else {
    if (is.null(object)) {
      stop("Provide `embedding` or a Seurat `object` + `reduction`.", call. = FALSE)
    }
    extract_reduction(object = object, reduction = reduction)
  }
  if (ncol(emb) < 2L) stop("`embedding` must have at least 2 columns.", call. = FALSE)
  if (is.null(rownames(emb))) {
    if (nrow(emb) != ncol(score$score)) {
      stop("Embedding has no rownames; nrow(embedding) must equal ncol(score).", call. = FALSE)
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
    value = as.numeric(score$score[pathway, ]),
    stringsAsFactors = FALSE
  )

  if (!is.null(split.by)) {
    grp <- resolve_meta_var(score$meta, split.by, "split.by")
    df$split <- as.factor(grp)
  }
  tp <- resolve_text_params(theme_params)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$dim1, y = ggplot2::.data$dim2, color = ggplot2::.data$value)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::labs(x = colnames(emb)[1], y = colnames(emb)[2], color = pathway) +
    do.call(gleam_theme, tp)

  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~ split)
  }
  p
}
