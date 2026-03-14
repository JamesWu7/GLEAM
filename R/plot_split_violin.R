#' Split violin plot for signature scores
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param pathway Legacy alias of `signature` (kept for backward compatibility).
#' @param x Grouping variable for x-axis.
#' @param split.by Variable used for split/fill.
#' @param palette Discrete palette name or custom colors.
#' @param alpha Violin transparency.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_split_violin <- function(
  score,
  signature = NULL,
  pathway = NULL,
  x,
  split.by,
  palette = "gleam_discrete",
  alpha = 0.75,
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature, pathway = pathway)

  xv <- resolve_meta_var(score$meta, x, "x")
  sv <- resolve_meta_var(score$meta, split.by, "split.by")
  x_chr <- as.character(xv)
  x_levels <- if (is.factor(xv)) levels(base::droplevels(xv)) else unique(x_chr)
  split_chr <- as.character(sv)
  split_levels <- if (is.factor(sv)) levels(base::droplevels(sv)) else unique(split_chr)
  if (!is.character(palette) && !is.null(names(palette)) && any(nzchar(names(palette)))) {
    pal_levels <- names(palette)[nzchar(names(palette))]
    split_levels <- c(pal_levels[pal_levels %in% split_chr], setdiff(split_levels, pal_levels))
  }
  df <- data.frame(
    x = factor(x_chr, levels = x_levels, ordered = TRUE),
    split = factor(split_chr, levels = split_levels, ordered = TRUE),
    value = as.numeric(score$score[signature, ]),
    stringsAsFactors = FALSE
  )
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$value, fill = .data$split)) +
    ggplot2::geom_violin(position = ggplot2::position_dodge(width = 0.85), alpha = alpha, trim = TRUE, color = "#1f2937", linewidth = 0.2) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.85)) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(title = paste("Split violin:", signature), x = "Group", y = "Signature score") +
    do.call(gleam_theme, tp)
}
