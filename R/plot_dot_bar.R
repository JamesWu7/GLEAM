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

  if (!is.null(pathway)) {
    mean_df <- mean_df[mean_df$signature %in% pathway, , drop = FALSE]
    frac_df <- frac_df[frac_df$signature %in% pathway, , drop = FALSE]
  }
  if (nrow(mean_df) < 1L) {
    stop("No signatures remain after filtering.", call. = FALSE)
  }

  key <- paste(mean_df$signature, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$signature, frac_df$group_key))
  mean_df$fraction <- frac_map[key]
  tp <- resolve_text_params(theme_params)

  cell_levels <- names(sort(tapply(
    X = mean_df$value,
    INDEX = as.character(mean_df[[celltype_col]]),
    FUN = function(x) mean(x, na.rm = TRUE)
  ), decreasing = TRUE))
  mean_df$celltype_plot <- factor(as.character(mean_df[[celltype_col]]), levels = rev(cell_levels))

  n_sig <- length(unique(as.character(mean_df$signature)))
  ncol_facet <- if (n_sig <= 2L) n_sig else min(3L, ceiling(sqrt(n_sig)))

  dot <- ggplot2::ggplot(mean_df, ggplot2::aes(
    x = .data$value,
    y = .data$celltype_plot,
    color = .data[[group_col]],
    size = .data$fraction
  )) +
    ggplot2::geom_point(alpha = 0.9, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    scale_gleam_color(color_palette, continuous = FALSE) +
    ggplot2::facet_wrap(~signature, scales = "free_x", ncol = ncol_facet) +
    ggplot2::labs(
      title = "Signature dot + bar summary",
      x = "Signature score",
      y = "Cell type",
      color = "Group",
      size = "Active fraction"
    ) +
    do.call(gleam_theme, tp)

  group_vals_raw <- mean_df[[group_col]]
  g_levels <- if (is.factor(group_vals_raw)) {
    levels(stats::droplevels(group_vals_raw))
  } else {
    unique(as.character(group_vals_raw))
  }
  g_levels <- g_levels[!is.na(g_levels) & nzchar(g_levels)]

  if (length(g_levels) < 2L) {
    return(dot)
  }
  if (length(g_levels) > 2L) {
    warning(
      "More than two groups detected; showing difference using the first two groups: ",
      g_levels[[1]], " and ", g_levels[[2]], ".",
      call. = FALSE
    )
  }
  g1 <- g_levels[[1]]
  g2 <- g_levels[[2]]

  sub <- mean_df[as.character(mean_df[[group_col]]) %in% c(g1, g2), c("signature", celltype_col, group_col, "value"), drop = FALSE]
  colnames(sub) <- c("signature", "celltype", "group", "value")
  agg <- stats::aggregate(
    x = sub$value,
    by = list(signature = sub$signature, celltype = sub$celltype, group = sub$group),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  colnames(agg)[colnames(agg) == "x"] <- "value"
  wide <- stats::reshape(
    agg,
    idvar = c("signature", "celltype"),
    timevar = "group",
    direction = "wide"
  )
  c1 <- paste0("value.", g1)
  c2 <- paste0("value.", g2)
  if (!all(c(c1, c2) %in% colnames(wide))) {
    return(dot)
  }
  wide$delta <- as.numeric(wide[[c2]]) - as.numeric(wide[[c1]])
  bar_df <- wide[, c("signature", "celltype", "delta"), drop = FALSE]
  bar_df <- bar_df[stats::complete.cases(bar_df), , drop = FALSE]
  if (is.null(bar_df) || nrow(bar_df) < 1L) {
    return(dot)
  }
  bar_df$signature <- factor(as.character(bar_df$signature), levels = unique(as.character(mean_df$signature)))
  bar_df$celltype_plot <- factor(as.character(bar_df$celltype), levels = levels(mean_df$celltype_plot))

  bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = .data$delta, y = .data$celltype_plot)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$delta >= 0), width = 0.65, color = NA) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "#2B83BA"), guide = "none") +
    ggplot2::facet_wrap(~signature, scales = "free_x", ncol = ncol_facet) +
    ggplot2::labs(
      title = paste0("Group difference (", g2, " - ", g1, ")"),
      x = "Delta signature score",
      y = "Cell type"
    ) +
    do.call(gleam_theme, tp) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(6, 6, 6, 0),
      legend.position = "none"
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' not installed; returning dot panel only.", call. = FALSE)
    return(dot)
  }

  dot + bar + patchwork::plot_layout(guides = "collect", widths = c(7, 4))
}
