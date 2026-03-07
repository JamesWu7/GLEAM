#' Heatmap of aggregated pathway scores
#'
#' @param score `gleam_score` object.
#' @param by Grouping metadata columns.
#' @param fun Aggregation function.
#' @param threshold Threshold used by `fraction`.
#' @param top_n Optional top pathways by variance.
#' @param palette Continuous palette name or custom colors.
#'
#' @return A `ggplot` object.
#' @export
plot_heatmap <- function(score, by, fun = c("mean", "median", "fraction"), threshold = 0, top_n = NULL, palette = "gleam_continuous") {
  fun <- match.arg(fun)
  agg <- aggregate_pathway(score, by = by, fun = fun, threshold = threshold, long = TRUE)

  if (!is.null(top_n) && is.numeric(top_n) && top_n > 0) {
    vv <- tapply(agg$value, agg$pathway, stats::var, na.rm = TRUE)
    keep <- names(sort(vv, decreasing = TRUE))[seq_len(min(top_n, length(vv)))]
    agg <- agg[agg$pathway %in% keep, , drop = FALSE]
  }

  agg$pathway <- factor(agg$pathway, levels = rev(unique(agg$pathway)))

  ggplot2::ggplot(agg, ggplot2::aes(x = group_key, y = pathway, fill = value)) +
    ggplot2::geom_tile() +
    scale_gleam_fill(palette = palette, continuous = TRUE) +
    ggplot2::labs(x = paste(by, collapse = ":"), y = "Pathway", title = "Pathway heatmap") +
    .theme_gleam()
}
