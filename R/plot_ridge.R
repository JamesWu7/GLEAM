#' Ridge plot for signature score distributions
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param group Group variable.
#' @param palette Discrete palette name or custom colors.
#' @param alpha Transparency of ridges.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_ridge <- function(score, signature = NULL, group, palette = "gleam_discrete", alpha = 0.75, theme_params = list()) {
  check_score_object(score)
  check_plot_dependency("plot_ridge", "ggridges")
  signature <- resolve_signature_arg(score, signature = signature)

  g <- resolve_meta_var(score$meta, group, "group")
  g_chr <- as.character(g)
  g_levels <- resolve_discrete_levels(g, palette = palette)
  pal_values <- resolve_discrete_palette_values(g_levels, palette = palette)
  df <- data.frame(
    group = factor(g_chr, levels = g_levels, ordered = TRUE),
    value = as.numeric(score$score[signature, ]),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$group, fill = .data$group)) +
    ggridges::geom_density_ridges(alpha = alpha, scale = 1.2, color = "#1f2937", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = pal_values, drop = FALSE) +
    ggplot2::labs(title = paste("Ridge plot:", signature), x = "Signature score", y = "Group") +
    do.call(gleam_theme, tp) +
    ggplot2::theme(legend.position = "none")
}
