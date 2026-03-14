#' Plot signature score over pseudotime
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param pathway Legacy alias of `signature` (kept for backward compatibility).
#' @param pseudotime Pseudotime source.
#' @param lineage Optional lineage source for coloring.
#' @param smooth Add smoothing line.
#' @param point_size Point size.
#' @param alpha Point alpha.
#' @param palette Discrete palette for lineages.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_pseudotime_score <- function(
  score,
  signature = NULL,
  pathway = NULL,
  pseudotime = NULL,
  lineage = NULL,
  smooth = TRUE,
  point_size = 1.1,
  alpha = 0.55,
  palette = "gleam_discrete",
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature, pathway = pathway)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  df <- data.frame(
    pseudotime = pt,
    score = as.numeric(score$score[signature, ]),
    lineage = as.factor(ln),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)
  n_lineage <- length(unique(as.character(df$lineage)))
  if (n_lineage < 1L) n_lineage <- 1L

  # Build exact-length discrete colors (Seurat-like behavior: enough colors for all groups).
  lineage_cols <- if (is.character(palette) && length(palette) == 1L) {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      base <- RColorBrewer::brewer.pal(9, "Set1")
      if (n_lineage <= length(base)) base[seq_len(n_lineage)] else grDevices::colorRampPalette(base)(n_lineage)
    } else {
      get_palette(name = palette, n = n_lineage, continuous = FALSE)
    }
  } else {
    pal <- as.character(palette)
    if (length(pal) >= n_lineage) pal[seq_len(n_lineage)] else grDevices::colorRampPalette(pal)(n_lineage)
  }
  names(lineage_cols) <- levels(df$lineage)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pseudotime, y = .data$score, color = .data$lineage)) +
    ggplot2::geom_point(alpha = alpha, size = point_size) +
    ggplot2::scale_color_manual(values = lineage_cols) +
    ggplot2::labs(title = paste("Pseudotime signature:", signature), x = "Pseudotime", y = "Signature score") +
    do.call(gleam_theme, tp)

  if (smooth) {
    p <- p + ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x)
  }
  p
}
