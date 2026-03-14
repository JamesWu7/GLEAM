#' Plot celltype comparison result
#'
#' @param x `gleam_test` object or result table.
#' @param top_n Number of signatures to display.
#' @param point_size Point size.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_celltype_compare <- function(x, top_n = 20, point_size = 2.5, theme_params = list()) {
  tbl <- if (inherits(x, "gleam_test")) x$table else x
  if (!is.data.frame(tbl)) stop("`x` must be a gleam_test object or data.frame.", call. = FALSE)

  ord <- order(abs(tbl$effect_size), decreasing = TRUE)
  sel <- tbl[utils::head(ord, min(top_n, nrow(tbl))), , drop = FALSE]
  sel$pathway <- factor(sel$pathway, levels = rev(sel$pathway))
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(sel, ggplot2::aes(x = .data$pathway, y = .data$effect_size, color = .data$p_adj)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Signature", y = "Effect size", title = "Celltype signature comparison") +
    do.call(gleam_theme, tp)
}
