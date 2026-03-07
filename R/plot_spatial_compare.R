#' Plot spatial comparison result
#'
#' @param x `gleam_test` object or result table.
#' @param top_n Number of pathways to display.
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_compare <- function(x, top_n = 20) {
  tbl <- if (inherits(x, "gleam_test") || inherits(x, "scpathway_test")) x$table else x
  if (!is.data.frame(tbl)) stop("`x` must be a gleam_test object or data.frame.", call. = FALSE)

  ord <- order(tbl$p_adj, decreasing = FALSE, na.last = TRUE)
  sel <- tbl[head(ord, min(top_n, nrow(tbl))), , drop = FALSE]
  sel$pathway <- factor(sel$pathway, levels = rev(sel$pathway))

  ggplot2::ggplot(sel, ggplot2::aes(x = pathway, y = effect_size, fill = direction)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Spatial differential pathways", x = "Pathway", y = "Effect size") +
    .theme_gleam()
}
