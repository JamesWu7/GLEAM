#' Plot signature score on embedding coordinates
#'
#' @param score `gleam_score` object.
#' @param pathway Signature name (legacy argument name).
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
  if (!pathway %in% rownames(score$score)) stop("`signature` not found in score matrix.", call. = FALSE)
  tp <- resolve_text_params(theme_params)

  # Seurat FeaturePlot-style rendering when Seurat object is supplied.
  if (is.null(embedding) && !is.null(object) && is_seurat_object(object) && requireNamespace("Seurat", quietly = TRUE)) {
    sig_col <- paste0("GLEAM_signature_", gsub("[^A-Za-z0-9_]+", "_", pathway))
    cells <- intersect(colnames(object), colnames(score$score))
    if (length(cells) < 1L) {
      stop("No overlapping cells between `object` and score matrix.", call. = FALSE)
    }
    vals <- rep(NA_real_, length(colnames(object)))
    names(vals) <- colnames(object)
    vals[cells] <- as.numeric(score$score[pathway, cells])
    object[[sig_col]] <- vals

    cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 128, continuous = TRUE) else palette
    p <- Seurat::FeaturePlot(
      object = object,
      features = sig_col,
      reduction = reduction,
      split.by = split.by,
      pt.size = point_size,
      cols = cols,
      order = TRUE
    )
    if ("patchwork" %in% class(p) && requireNamespace("patchwork", quietly = TRUE)) {
      p <- p & do.call(gleam_theme, tp)
    } else if (inherits(p, "ggplot")) {
      p <- p + do.call(gleam_theme, tp)
    }
    return(p)
  }

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
  # draw high-value cells on top to mimic FeaturePlot readability
  df <- df[order(df$value), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$value)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::labs(
      title = paste("Signature:", pathway),
      x = colnames(emb)[1],
      y = colnames(emb)[2],
      color = "Signature score"
    ) +
    do.call(gleam_theme, tp)

  if (!is.null(split.by)) {
    p <- p + ggplot2::facet_wrap(~ split)
  }
  if (tolower(reduction) %in% c("umap", "tsne")) {
    p <- p + ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
  }
  p
}
