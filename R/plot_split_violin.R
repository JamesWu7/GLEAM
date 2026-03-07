#' Split violin plot for pathway scores
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param x Grouping variable for x-axis.
#' @param split.by Variable used for split/fill.
#'
#' @return A `ggplot` object.
#' @export
plot_split_violin <- function(score, pathway, x, split.by) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  xv <- resolve_meta_var(score$meta, x, "x")
  sv <- resolve_meta_var(score$meta, split.by, "split.by")
  df <- data.frame(x = as.factor(xv), split = as.factor(sv), value = as.numeric(score$score[pathway, ]), stringsAsFactors = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = value, fill = split)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.85), alpha = 0.75, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.85)) +
    ggplot2::labs(title = paste("Split violin:", pathway), x = as.character(substitute(x)), y = "Pathway score") +
    .theme_gleam()
}
