#' Dot plot of pathway summaries
#'
#' Dot size represents fraction of cells above threshold,
#' and dot color represents mean pathway score.
#'
#' @param score `gleam_score` object.
#' @param by Grouping metadata columns.
#' @param threshold Threshold for active fraction.
#' @param palette Continuous palette name or custom colors.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_dot <- function(score, by, threshold = 0, palette = "gleam_continuous", theme_params = list()) {
  mean_df <- aggregate_pathway(score, by = by, fun = "mean", long = TRUE)
  frac_df <- aggregate_pathway(score, by = by, fun = "fraction", threshold = threshold, long = TRUE)

  key <- paste(mean_df$pathway, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$pathway, frac_df$group_key))
  mean_df$fraction <- frac_map[key]
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(mean_df, ggplot2::aes(x = .data$group_key, y = .data$pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$fraction, color = .data$value), alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::labs(x = paste(by, collapse = ":"), y = "Pathway", title = "Pathway dot plot") +
    do.call(gleam_theme, tp)
}
