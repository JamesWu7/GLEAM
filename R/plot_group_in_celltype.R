#' Plot group comparison within a fixed celltype
#'
#' @param x `gleam_test` object or result table.
#' @param top_n Number of pathways to display.
#'
#' @return A `ggplot` object.
#' @export
plot_group_in_celltype <- function(x, top_n = 20) {
  tbl <- if (inherits(x, "gleam_test") || inherits(x, "scpathway_test")) x$table else x
  if (!is.data.frame(tbl)) stop("`x` must be a gleam_test object or data.frame.", call. = FALSE)

  ord <- order(tbl$p_adj, decreasing = FALSE, na.last = TRUE)
  sel <- tbl[utils::head(ord, min(top_n, nrow(tbl))), , drop = FALSE]
  sel$pathway <- factor(sel$pathway, levels = rev(sel$pathway))

  ggplot2::ggplot(sel, ggplot2::aes(x = ggplot2::.data$pathway, y = ggplot2::.data$effect_size, fill = ggplot2::.data$direction)) +
    ggplot2::geom_col(alpha = 0.85, width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Pathway", y = "Effect size", title = "Group difference within celltype") +
    .theme_gleam()
}
