#' Dot + bar summary plot for signature groups
#'
#' @param score `gleam_score` object.
#' @param by Grouping metadata columns.
#' @param threshold Threshold for active fraction.
#' @param pathway Optional signature filter (legacy argument name).
#' @param color_palette Continuous palette for mean score.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_dot_bar <- function(score, by, threshold = 0, pathway = NULL, color_palette = "gleam_continuous", theme_params = list()) {
  mean_df <- aggregate_signature(score, by = by, fun = "mean", long = TRUE)
  frac_df <- aggregate_signature(score, by = by, fun = "fraction", threshold = threshold, long = TRUE)
  mean_df$signature <- mean_df$pathway
  frac_df$signature <- frac_df$pathway

  if (!is.null(pathway)) {
    mean_df <- mean_df[mean_df$signature %in% pathway, , drop = FALSE]
    frac_df <- frac_df[frac_df$signature %in% pathway, , drop = FALSE]
  }

  key <- paste(mean_df$signature, mean_df$group_key)
  frac_map <- stats::setNames(frac_df$value, paste(frac_df$signature, frac_df$group_key))
  mean_df$fraction <- frac_map[key]
  tp <- resolve_text_params(theme_params)

  dot <- ggplot2::ggplot(mean_df, ggplot2::aes(x = .data$group_key, y = .data$signature)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$fraction, color = .data$value), alpha = 0.9) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    scale_gleam_color(color_palette, continuous = TRUE) +
    ggplot2::labs(
      title = "Signature dot + bar summary",
      x = paste(by, collapse = ":"),
      y = "Signature",
      color = "Mean signature score",
      size = "Active fraction"
    ) +
    do.call(gleam_theme, tp)

  # Build side bar panel when at least two groups are available.
  group_col <- if ("group" %in% colnames(mean_df)) "group" else by[[1]]
  bar_df <- NULL
  if (!is.null(group_col) && group_col %in% colnames(mean_df)) {
    g_levels <- unique(as.character(mean_df[[group_col]]))
    if (length(g_levels) >= 2L) {
      agg <- stats::aggregate(
        x = mean_df$value,
        by = list(
          signature = mean_df$signature,
          group = as.character(mean_df[[group_col]])
        ),
        FUN = function(x) mean(x, na.rm = TRUE)
      )
      colnames(agg)[colnames(agg) == "x"] <- "value"
      wide <- stats::reshape(
        agg,
        idvar = "signature",
        timevar = "group",
        direction = "wide"
      )
      g1 <- paste0("value.", g_levels[[1]])
      g2 <- paste0("value.", g_levels[[2]])
      if (all(c(g1, g2) %in% colnames(wide))) {
        denom <- pmax(abs(wide[[g1]]), 1e-8)
        wide$percent_increase <- (wide[[g2]] - wide[[g1]]) / denom * 100
        bar_df <- wide[, c("signature", "percent_increase"), drop = FALSE]
      }
    }
  }

  if (is.null(bar_df) || nrow(bar_df) < 1L) {
    return(dot)
  }

  bar_df$signature <- factor(as.character(bar_df$signature), levels = rev(unique(as.character(mean_df$signature))))
  bar <- ggplot2::ggplot(bar_df, ggplot2::aes(x = .data$percent_increase, y = .data$signature)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$percent_increase >= 0), width = 0.65, color = NA) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D73027", `FALSE` = "#2B83BA"), guide = "none") +
    ggplot2::labs(title = "% Increase", x = NULL, y = NULL) +
    do.call(gleam_theme, tp) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 6, 6, 0)
    )

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' not installed; returning dot panel only.", call. = FALSE)
    return(dot)
  }

  dot + bar + patchwork::plot_layout(guides = "collect", widths = c(7, 3))
}
