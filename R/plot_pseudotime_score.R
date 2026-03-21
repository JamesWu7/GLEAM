#' Plot signature score over pseudotime
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
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
  pseudotime = NULL,
  lineage = NULL,
  smooth = TRUE,
  point_size = 1.1,
  alpha = 0.55,
  palette = "gleam_discrete",
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  df <- data.frame(
    pseudotime = pt,
    score = as.numeric(score$score[signature, ]),
    lineage = as.factor(ln),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)
  lineage_levels <- resolve_discrete_levels(df$lineage, palette = palette)
  lineage_cols <- resolve_discrete_palette_values(lineage_levels, palette = palette)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pseudotime, y = .data$score, color = .data$lineage)) +
    ggplot2::geom_point(alpha = alpha, size = point_size) +
    ggplot2::scale_color_manual(values = lineage_cols, drop = FALSE) +
    ggplot2::labs(title = paste("Pseudotime signature:", signature), x = "Pseudotime", y = "Signature score") +
    do.call(gleam_theme, tp)

  if (smooth) {
    p <- p + ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x)
  }
  p
}
