#' Plot pathway score over pseudotime
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
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
  pathway,
  pseudotime = NULL,
  lineage = NULL,
  smooth = TRUE,
  point_size = 1.1,
  alpha = 0.55,
  palette = "gleam_discrete",
  theme_params = list()
) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  df <- data.frame(
    pseudotime = pt,
    score = as.numeric(score$score[pathway, ]),
    lineage = as.factor(ln),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$pseudotime, y = ggplot2::.data$score, color = ggplot2::.data$lineage)) +
    ggplot2::geom_point(alpha = alpha, size = point_size) +
    scale_gleam_color(palette = palette, continuous = FALSE) +
    ggplot2::labs(title = paste("Pseudotime score:", pathway), x = "Pseudotime", y = "Pathway score") +
    do.call(gleam_theme, tp)

  if (smooth) {
    p <- p + ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x)
  }
  p
}
