#' Split violin plot for pathway scores
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param x Grouping variable for x-axis.
#' @param split.by Variable used for split/fill.
#' @param palette Discrete palette name or custom colors.
#' @param alpha Violin transparency.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_split_violin <- function(
  score,
  pathway,
  x,
  split.by,
  palette = "gleam_discrete",
  alpha = 0.75,
  theme_params = list()
) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  xv <- resolve_meta_var(score$meta, x, "x")
  sv <- resolve_meta_var(score$meta, split.by, "split.by")
  df <- data.frame(x = as.factor(xv), split = as.factor(sv), value = as.numeric(score$score[pathway, ]), stringsAsFactors = FALSE)
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$value, fill = .data$split)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.85), alpha = alpha, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.85)) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(title = paste("Split violin:", pathway), x = as.character(substitute(x)), y = "Pathway score") +
    do.call(gleam_theme, tp)
}
