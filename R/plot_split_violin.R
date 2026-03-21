#' Split violin plot for signature scores
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
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
  signature = NULL,
  x,
  split.by,
  palette = "gleam_discrete",
  alpha = 0.75,
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature)

  xv <- resolve_meta_var(score$meta, x, "x")
  sv <- resolve_meta_var(score$meta, split.by, "split.by")
  x_chr <- as.character(xv)
  x_levels <- resolve_discrete_levels(xv, palette = "gleam_discrete")
  split_chr <- as.character(sv)
  split_levels <- resolve_discrete_levels(sv, palette = palette)
  pal_values <- resolve_discrete_palette_values(split_levels, palette = palette)
  df <- data.frame(
    x = factor(x_chr, levels = x_levels, ordered = TRUE),
    split = factor(split_chr, levels = split_levels, ordered = TRUE),
    value = as.numeric(score$score[signature, ]),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$value, fill = .data$split)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.85), alpha = alpha, trim = TRUE, color = "#1f2937", linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.85)) +
    ggplot2::scale_fill_manual(values = pal_values, drop = FALSE) +
    ggplot2::labs(title = paste("Split violin:", signature), x = "Group", y = "Signature score") +
    do.call(gleam_theme, tp)
}
