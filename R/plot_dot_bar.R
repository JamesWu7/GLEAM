#' Dot + bar summary plot for pathway groups
#'
#' @param score `gleam_score` object.
#' @param by Grouping metadata columns.
#' @param threshold Threshold for active fraction.
#' @param pathway Optional pathway filter.
#'
#' @return A `ggplot` object.
#' @export
plot_dot_bar <- function(score, by, threshold = 0, pathway = NULL) {
  mean_df <- aggregate_pathway(score, by = by, fun = "mean", long = TRUE)
  frac_df <- aggregate_pathway(score, by = by, fun = "fraction", threshold = threshold, long = TRUE)

  if (!is.null(pathway)) {
    mean_df <- mean_df[mean_df$pathway %in% pathway, , drop = FALSE]
    frac_df <- frac_df[frac_df$pathway %in% pathway, , drop = FALSE]
  }

  key <- paste(mean_df$pathway, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$pathway, frac_df$group_key))
  mean_df$fraction <- frac_map[key]

  ggplot2::ggplot(mean_df, ggplot2::aes(x = ggplot2::.data$group_key, y = ggplot2::.data$pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = ggplot2::.data$fraction, color = ggplot2::.data$value), alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    scale_gleam_color("gleam_continuous", continuous = TRUE) +
    ggplot2::labs(title = "Dot-bar summary", x = paste(by, collapse = ":"), y = "Pathway") +
    .theme_gleam()
}
