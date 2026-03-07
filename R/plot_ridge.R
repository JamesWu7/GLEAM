#' Ridge plot for pathway score distributions
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param group Group variable.
#' @param palette Discrete palette name or custom colors.
#' @param alpha Transparency of ridges.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_ridge <- function(score, pathway, group, palette = "gleam_discrete", alpha = 0.75, theme_params = list()) {
  check_score_object(score)
  check_plot_dependency("plot_ridge", "ggridges")
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  g <- resolve_meta_var(score$meta, group, "group")
  df <- data.frame(group = as.factor(g), value = as.numeric(score$score[pathway, ]), stringsAsFactors = FALSE)
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$value, y = ggplot2::.data$group, fill = ggplot2::.data$group)) +
    ggridges::geom_density_ridges(alpha = alpha, scale = 1.2, color = "#1f2937", linewidth = 0.25) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(title = paste("Ridge plot:", pathway), x = "Pathway score", y = "Group") +
    do.call(gleam_theme, tp)
}
