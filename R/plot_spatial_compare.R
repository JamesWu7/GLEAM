#' Plot spatial comparison result
#'
#' @param x `gleam_test` object or result table.
#' @param top_n Number of pathways to display.
#' @param alpha Bar alpha.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_compare <- function(x, top_n = 20, alpha = 0.85, theme_params = list()) {
  tbl <- if (inherits(x, "gleam_test") || inherits(x, "scpathway_test")) x$table else x
  if (!is.data.frame(tbl)) stop("`x` must be a gleam_test object or data.frame.", call. = FALSE)

  ord <- order(tbl$p_adj, decreasing = FALSE, na.last = TRUE)
  sel <- tbl[utils::head(ord, min(top_n, nrow(tbl))), , drop = FALSE]
  sel$pathway <- factor(sel$pathway, levels = rev(sel$pathway))
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(sel, ggplot2::aes(x = .data$pathway, y = .data$effect_size, fill = .data$direction)) +
    ggplot2::geom_col(alpha = alpha) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Spatial differential pathways", x = "Pathway", y = "Effect size") +
    do.call(gleam_theme, tp)
}
