#' Violin plot of signature scores
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param pathway Legacy alias of `signature` (kept for backward compatibility).
#' @param group Group variable (metadata column name or vector).
#' @param celltype Optional celltype filter (single label or vector).
#' @param trim Passed to `geom_violin()`.
#' @param palette Discrete palette name or custom colors.
#' @param point_size Point size for boxplot outlier overlay.
#' @param alpha Violin transparency.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_violin <- function(
  score,
  signature = NULL,
  pathway = NULL,
  group,
  celltype = NULL,
  trim = TRUE,
  palette = "gleam_discrete",
  point_size = 1.4,
  alpha = 0.7,
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature, pathway = pathway)

  meta <- score$meta
  g <- resolve_meta_var(meta, group, "group")
  y <- as.numeric(score$score[signature, ])

  keep <- rep(TRUE, length(y))
  if (!is.null(celltype)) {
    if (length(celltype) == 1L && is.character(celltype) && celltype %in% colnames(meta)) {
      keep <- !is.na(as.character(meta[[celltype]]))
    } else if (length(celltype) == 1L) {
      if (!"celltype" %in% colnames(meta)) {
        stop("`celltype` filter supplied but metadata has no 'celltype' column.", call. = FALSE)
      }
      keep <- as.character(meta$celltype) %in% as.character(celltype)
    } else {
      if (!"celltype" %in% colnames(meta)) {
        stop("`celltype` filter supplied but metadata has no 'celltype' column.", call. = FALSE)
      }
      keep <- as.character(meta$celltype) %in% as.character(celltype)
    }
  }

  g_keep <- g[keep]
  g_chr <- as.character(g_keep)
  g_levels <- if (is.factor(g_keep)) levels(base::droplevels(g_keep)) else unique(g_chr)
  if (!is.character(palette) && !is.null(names(palette)) && any(nzchar(names(palette)))) {
    pal_levels <- names(palette)[nzchar(names(palette))]
    g_levels <- c(pal_levels[pal_levels %in% g_chr], setdiff(g_levels, pal_levels))
  }
  df <- data.frame(
    group = factor(g_chr, levels = g_levels, ordered = TRUE),
    score = y[keep],
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$score, fill = .data$group)) +
    ggplot2::geom_violin(trim = trim, alpha = alpha, color = NA) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8, linewidth = 0.3) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(x = "Group", y = "Signature score", title = paste("Signature score:", signature)) +
    do.call(gleam_theme, tp)

  if (isTRUE(point_size > 0)) {
    p <- p + ggplot2::geom_jitter(
      width = 0.08,
      size = point_size,
      alpha = min(1, alpha + 0.1),
      color = "#111827"
    )
  }
  p
}
