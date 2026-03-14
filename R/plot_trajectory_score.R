#' Plot signature score on trajectory embedding
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param pathway Legacy alias of `signature` (kept for backward compatibility).
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
  signature = NULL,
  pathway = NULL,
  embeddings = NULL,
  reduction = "umap",
  object = NULL,
  point_size = 1.1,
  alpha = 0.9,
  palette = "gleam_continuous",
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature, pathway = pathway)

  if (is.null(embeddings) && is.null(object)) {
    stop("Provide `embeddings` or a Seurat `object`.", call. = FALSE)
  }

  emb <- if (!is.null(embeddings)) {
    as.matrix(embeddings)
  } else {
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
    cell_id = colnames(score$score),
    dim1 = emb[, 1],
    dim2 = emb[, 2],
    score = as.numeric(score$score[signature, ]),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  line_df <- NULL
  end_df <- NULL
  meta <- score$meta
  if (!is.null(meta) && all(c("pseudotime", "lineage") %in% colnames(meta))) {
    if (is.null(rownames(meta)) && nrow(meta) == ncol(score$score)) {
      rownames(meta) <- colnames(score$score)
    }
    if (!is.null(rownames(meta)) && all(colnames(score$score) %in% rownames(meta))) {
      meta <- meta[colnames(score$score), , drop = FALSE]
      traj_df <- data.frame(
        cell_id = rownames(meta),
        dim1 = df$dim1,
        dim2 = df$dim2,
        pseudotime = suppressWarnings(as.numeric(meta$pseudotime)),
        lineage = as.character(meta$lineage),
        stringsAsFactors = FALSE
      )
      traj_df <- traj_df[is.finite(traj_df$pseudotime) & !is.na(traj_df$lineage) & nzchar(traj_df$lineage), , drop = FALSE]
      if (nrow(traj_df) > 0L) {
        lines_list <- lapply(split(traj_df, traj_df$lineage), function(d) {
          if (nrow(d) < 20L) return(NULL)
          d <- d[order(d$pseudotime), , drop = FALSE]
          n_bins <- min(35L, max(8L, floor(nrow(d) / 35L)))
          brks <- unique(as.numeric(stats::quantile(
            d$pseudotime,
            probs = seq(0, 1, length.out = n_bins + 1L),
            na.rm = TRUE,
            names = FALSE,
            type = 7
          )))
          if (length(brks) < 3L) return(NULL)
          d$bin <- cut(d$pseudotime, breaks = brks, include.lowest = TRUE, labels = FALSE)
          agg <- stats::aggregate(
            cbind(dim1, dim2, pseudotime) ~ bin,
            data = d,
            FUN = function(x) stats::median(x, na.rm = TRUE)
          )
          agg <- agg[order(agg$pseudotime), , drop = FALSE]
          if (nrow(agg) < 3L) return(NULL)
          agg$lineage <- as.character(d$lineage[[1]])
          agg$step <- seq_len(nrow(agg))
          agg
        })
        lines_list <- Filter(Negate(is.null), lines_list)
        if (length(lines_list) > 0L) {
          line_df <- do.call(rbind, lines_list)
          end_df <- do.call(rbind, lapply(split(line_df, line_df$lineage), function(d) {
            d <- d[order(d$step), , drop = FALSE]
            d[nrow(d), c("lineage", "dim1", "dim2"), drop = FALSE]
          }))
        }
      }
    }
  }

  df <- df[order(df$score), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$dim1, y = .data$dim2, color = .data$score)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::labs(title = paste("Trajectory signature map:", signature), x = colnames(emb)[1], y = colnames(emb)[2], color = "Signature score") +
    do.call(gleam_theme, tp)

  # Overlay lineage curves in a slingshot-like style when pseudotime+lineage exist.
  if (!is.null(line_df) && nrow(line_df) > 0L) {
    p <- p +
      ggplot2::geom_path(
        data = line_df,
        ggplot2::aes(x = .data$dim1, y = .data$dim2, group = .data$lineage),
        inherit.aes = FALSE,
        color = "#111827",
        linewidth = 1.15,
        alpha = 0.92,
        lineend = "round"
      ) +
      ggplot2::geom_point(
        data = end_df,
        ggplot2::aes(x = .data$dim1, y = .data$dim2),
        inherit.aes = FALSE,
        shape = 21,
        size = 2.8,
        stroke = 0.35,
        fill = "#111827",
        color = "white"
      ) +
      ggplot2::geom_text(
        data = end_df,
        ggplot2::aes(x = .data$dim1, y = .data$dim2, label = .data$lineage),
        inherit.aes = FALSE,
        color = "#111827",
        size = 3.1,
        hjust = -0.12,
        vjust = -0.58
      )
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
