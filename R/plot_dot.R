#' Dot plot of pathway summaries
#'
#' Dot size represents fraction of cells above threshold,
#' and dot color represents mean pathway score.
#'
#' @param score `scpathway_score` object.
#' @param by Grouping metadata columns.
#' @param threshold Threshold for active fraction.
#'
#' @return A `ggplot` object.
#' @export
plot_dot <- function(score, by, threshold = 0) {
  mean_df <- aggregate_pathway(score, by = by, fun = "mean", long = TRUE)
  frac_df <- aggregate_pathway(score, by = by, fun = "fraction", threshold = threshold, long = TRUE)

  key <- paste(mean_df$pathway, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$pathway, frac_df$group_key))
  mean_df$fraction <- frac_map[key]

  ggplot2::ggplot(mean_df, ggplot2::aes(x = group_key, y = pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = fraction, color = value), alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    ggplot2::scale_color_gradient2(low = "#3b4cc0", mid = "#f7f7f7", high = "#b40426", midpoint = 0) +
    ggplot2::labs(x = paste(by, collapse = ":"), y = "Pathway", title = "Pathway dot plot") +
    .theme_scpathway()
}
