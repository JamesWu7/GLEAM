#' Plot pathway score on embedding coordinates
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param embedding Embedding matrix with at least 2 columns.
#' @param split.by Optional metadata column for faceting.
#'
#' @return A `ggplot` object.
#' @export
plot_embedding_score <- function(score, pathway, embedding, split.by = NULL) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)
  emb <- as.matrix(embedding)
  if (ncol(emb) < 2L) stop("`embedding` must have at least 2 columns.", call. = FALSE)

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

  p <- ggplot2::ggplot(df, ggplot2::aes(dim1, dim2, color = value)) +
    ggplot2::geom_point(size = 1.1, alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#3b4cc0", mid = "#f7f7f7", high = "#b40426", midpoint = 0) +
    ggplot2::labs(x = colnames(emb)[1], y = colnames(emb)[2], color = pathway) +
    .theme_gleam()

  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~ split)
  }
  p
}
