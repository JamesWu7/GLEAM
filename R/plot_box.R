#' Box plot of sample-level signature scores
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param group Group variable (metadata column name or vector).
#' @param sample Sample variable (metadata column name or vector).
#' @param palette Discrete palette name or custom colors.
#' @param point_size Jitter point size.
#' @param alpha Box alpha.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_box <- function(
  score,
  signature,
  group,
  sample,
  palette = "gleam_discrete",
  point_size = 1.8,
  alpha = 0.7,
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature)

  meta <- score$meta
  g <- resolve_meta_var(meta, group, "group")
  s <- resolve_meta_var(meta, sample, "sample")
  y <- as.numeric(score$score[signature, ])

  df <- data.frame(sample = as.character(s), group = as.character(g), score = y, stringsAsFactors = FALSE)
  samp <- stats::aggregate(score ~ sample + group, data = df, FUN = mean)
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(samp, ggplot2::aes(x = .data$group, y = .data$score, fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = alpha, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.12, size = point_size, alpha = min(1, alpha + 0.1)) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(x = "Group", y = "Signature score", title = paste("Sample-level signature:", signature)) +
    do.call(gleam_theme, tp)
}
