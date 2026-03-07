#' Plot celltype comparison result
#'
#' @param x `gleam_test` object or result table.
#' @param top_n Number of pathways to display.
#'
#' @return A `ggplot` object.
#' @export
plot_celltype_compare <- function(x, top_n = 20) {
  tbl <- if (inherits(x, "gleam_test") || inherits(x, "scpathway_test")) x$table else x
  if (!is.data.frame(tbl)) stop("`x` must be a gleam_test object or data.frame.", call. = FALSE)

  ord <- order(abs(tbl$effect_size), decreasing = TRUE)
  sel <- tbl[utils::head(ord, min(top_n, nrow(tbl))), , drop = FALSE]
  sel$pathway <- factor(sel$pathway, levels = rev(sel$pathway))

  ggplot2::ggplot(sel, ggplot2::aes(x = ggplot2::.data$pathway, y = ggplot2::.data$effect_size, color = ggplot2::.data$p_adj)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Pathway", y = "Effect size", title = "Celltype pathway comparison") +
    .theme_gleam()
}
