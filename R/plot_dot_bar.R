#' Dot + bar summary plot for signature groups
#'
#' @param score `gleam_score` object.
#' @param by Grouping metadata columns. The first element is interpreted as
#'   `group`, and the second as `celltype`.
#' @param threshold Threshold for active fraction.
#' @param pathway Optional signature filter (legacy argument name).
#' @param color_palette Continuous palette for mean score.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_dot_bar <- function(score, by, threshold = 0, pathway = NULL, color_palette = "gleam_continuous", theme_params = list()) {
  if (!is.character(by) || length(by) < 2L) {
    stop("`by` must be character columns in order: c(group_col, celltype_col).", call. = FALSE)
  }
  group_col <- by[[1]]
  celltype_col <- by[[2]]

  mean_df <- aggregate_signature(score, by = by, fun = "mean", long = TRUE)
  frac_df <- aggregate_signature(score, by = by, fun = "fraction", threshold = threshold, long = TRUE)
  if (!all(c(group_col, celltype_col) %in% colnames(mean_df))) {
    stop("`by` columns must exist in `score$meta` and be present in aggregated output.", call. = FALSE)
  }

  mean_df$signature <- mean_df$pathway
  frac_df$signature <- frac_df$pathway

  sig_use <- if (is.null(pathway)) rownames(score$score)[[1]] else as.character(pathway)[[1]]
  if (!is.null(pathway) && length(pathway) > 1L) {
    warning("`plot_dot_bar()` now supports one signature at a time; using the first signature only.", call. = FALSE)
  }

  mean_df <- mean_df[as.character(mean_df$signature) == sig_use, , drop = FALSE]
  frac_df <- frac_df[as.character(frac_df$signature) == sig_use, , drop = FALSE]
  if (nrow(mean_df) < 1L) {
    stop("Selected signature not found in aggregated output.", call. = FALSE)
  }

  key <- paste(mean_df$signature, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$signature, frac_df$group_key))
  mean_df$fraction <- as.numeric(frac_map[key])
  tp <- resolve_text_params(theme_params)

  cell_raw <- mean_df[[celltype_col]]
  if (is.factor(cell_raw)) {
    base_levels <- levels(base::droplevels(cell_raw))
  } else {
    base_levels <- names(sort(tapply(
      X = mean_df$value,
      INDEX = as.character(cell_raw),
      FUN = function(x) mean(x, na.rm = TRUE)
    ), decreasing = TRUE))
  }
  display_levels <- base_levels
  y_levels <- rev(display_levels)
  mean_df$celltype_plot <- factor(as.character(mean_df[[celltype_col]]), levels = display_levels, ordered = TRUE)

  dot_df <- stats::aggregate(
    x = list(value = mean_df$value, fraction = mean_df$fraction),
    by = list(celltype = mean_df$celltype_plot),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  dot_df$celltype <- factor(as.character(dot_df$celltype), levels = levels(mean_df$celltype_plot), ordered = TRUE)
  dot_df$signature <- sig_use

  dot <- ggplot2::ggplot(dot_df, ggplot2::aes(
    x = .data$signature,
    y = .data$celltype,
    color = .data$value,
    size = .data$fraction
  )) +
    ggplot2::geom_point(alpha = 0.95, stroke = 0.25) +
    scale_gleam_color(color_palette, continuous = TRUE) +
    ggplot2::scale_size_continuous(range = c(2.6, 10.8)) +
    ggplot2::scale_y_discrete(limits = y_levels) +
    ggplot2::labs(
      title = "Signature score and active fraction",
      x = "Signature",
      y = "Cell type",
      color = "Mean score",
      size = "Active fraction"
    ) +
    do.call(gleam_theme, tp) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 8, 6, 8),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1, vjust = 1),
      legend.box = "vertical",
      legend.spacing.y = ggplot2::unit(2, "mm")
    )

  group_vals_raw <- mean_df[[group_col]]
  g_levels <- if (is.factor(group_vals_raw)) levels(base::droplevels(group_vals_raw)) else unique(as.character(group_vals_raw))
  g_levels <- g_levels[!is.na(g_levels) & nzchar(g_levels)]
  if (length(g_levels) < 2L) {
    return(dot)
  }
  if (length(g_levels) > 2L) {
    warning(
      "More than two groups detected; right-panel delta uses the first two groups: ",
      g_levels[[1]], " and ", g_levels[[2]], ".",
      call. = FALSE
    )
  }
  g1 <- g_levels[[1]]
  g2 <- g_levels[[2]]

  sub <- mean_df[as.character(mean_df[[group_col]]) %in% c(g1, g2), c(celltype_col, group_col, "value"), drop = FALSE]
  colnames(sub) <- c("celltype", "group", "value")
  agg <- stats::aggregate(
    x = sub$value,
    by = list(celltype = sub$celltype, group = sub$group),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  colnames(agg)[colnames(agg) == "x"] <- "value"
  wide <- stats::reshape(agg, idvar = "celltype", timevar = "group", direction = "wide")
  c1 <- paste0("value.", g1)
  c2 <- paste0("value.", g2)
  if (!all(c(c1, c2) %in% colnames(wide))) {
    return(dot)
  }

  bar_df <- wide[, c("celltype", c1, c2), drop = FALSE]
  bar_df <- bar_df[stats::complete.cases(bar_df), , drop = FALSE]
  if (nrow(bar_df) < 1L) {
    return(dot)
  }
  bar_df$delta <- as.numeric(bar_df[[c2]]) - as.numeric(bar_df[[c1]])
  bar_df$celltype_plot <- factor(as.character(bar_df$celltype), levels = levels(mean_df$celltype_plot), ordered = TRUE)
  bar_df <- bar_df[order(bar_df$celltype_plot), , drop = FALSE]
  bar_df$direction <- factor(ifelse(bar_df$delta >= 0, "increase", "decrease"), levels = c("decrease", "increase"))
  ct_levels <- levels(bar_df$celltype_plot)
  ct_cols_all <- get_palette("gleam_discrete", n = length(base_levels), continuous = FALSE)
  names(ct_cols_all) <- base_levels
  ct_cols <- ct_cols_all[ct_levels]

  bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = .data$delta, y = .data$celltype_plot)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$celltype_plot), width = 0.64, color = NA, alpha = 0.9) +
    ggplot2::geom_path(ggplot2::aes(group = 1), linewidth = 0.55, color = "#111827", alpha = 0.9) +
    ggplot2::geom_point(size = 1.8, color = "#111827") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = ct_cols, drop = FALSE) +
    ggplot2::scale_y_discrete(limits = y_levels) +
    ggplot2::labs(
      title = paste0("Group difference curve (", g2, " - ", g1, ")"),
      x = "Delta signature score",
      y = "Cell type"
    ) +
    do.call(gleam_theme, tp) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 6, 6, 2),
      legend.position = "none"
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' not installed; returning dot panel only.", call. = FALSE)
    return(dot)
  }

  dot + bar + patchwork::plot_layout(guides = "collect", widths = c(3, 5))
}
