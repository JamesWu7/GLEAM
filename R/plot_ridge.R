#' Ridge plot for pathway score distributions
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param group Group variable.
#'
#' @return A `ggplot` object.
#' @export
plot_ridge <- function(score, pathway, group) {
  check_score_object(score)
  check_plot_dependency("plot_ridge", "ggridges")
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  g <- resolve_meta_var(score$meta, group, "group")
  df <- data.frame(group = as.factor(g), value = as.numeric(score$score[pathway, ]), stringsAsFactors = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$value, y = ggplot2::.data$group, fill = ggplot2::.data$group)) +
    ggridges::geom_density_ridges(alpha = 0.75, scale = 1.2) +
    ggplot2::labs(title = paste("Ridge plot:", pathway), x = "Pathway score", y = "Group") +
    .theme_gleam()
}
