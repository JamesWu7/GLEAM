#' Box plot of sample-level pathway scores
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param group Group variable (metadata column name or vector).
#' @param sample Sample variable (metadata column name or vector).
#' @param palette Discrete palette name or custom colors.
#'
#' @return A `ggplot` object.
#' @export
plot_box <- function(score, pathway, group, sample, palette = "gleam_discrete") {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) {
    stop("`pathway` not found in score matrix.", call. = FALSE)
  }

  meta <- score$meta
  g <- resolve_meta_var(meta, group, "group")
  s <- resolve_meta_var(meta, sample, "sample")
  y <- as.numeric(score$score[pathway, ])

  df <- data.frame(sample = as.character(s), group = as.character(g), score = y, stringsAsFactors = FALSE)
  samp <- stats::aggregate(score ~ sample + group, data = df, FUN = mean)

  ggplot2::ggplot(samp, ggplot2::aes(x = ggplot2::.data$group, y = ggplot2::.data$score, fill = ggplot2::.data$group)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = 1.8, alpha = 0.8) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(x = "Group", y = pathway, title = paste("Sample-level:", pathway)) +
    .theme_gleam()
}
