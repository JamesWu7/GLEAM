#' Pseudobulk boxplot for signature scores
#'
#' @param score `gleam_score` object.
#' @param signature Signature name.
#' @param group Group variable.
#' @param sample Sample variable.
#' @param celltype Optional celltype variable.
#' @param palette Discrete palette name or custom colors.
#' @param point_size Jitter point size.
#' @param alpha Box alpha.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_pseudobulk_box <- function(
  score,
  signature,
  group,
  sample,
  celltype = NULL,
  palette = "gleam_discrete",
  point_size = 1.6,
  alpha = 0.75,
  theme_params = list()
) {
  check_score_object(score)
  signature <- resolve_signature_arg(score, signature = signature)

  g <- resolve_meta_var(score$meta, group, "group")
  s <- resolve_meta_var(score$meta, sample, "sample")
  ctv <- if (!is.null(celltype)) resolve_meta_var(score$meta, celltype, "celltype") else rep("all", ncol(score$score))

  df <- data.frame(
    sample = as.character(s),
    group = as.character(g),
    celltype = as.character(ctv),
    value = as.numeric(score$score[signature, ]),
    stringsAsFactors = FALSE
  )

  pb <- stats::aggregate(value ~ sample + group + celltype, data = df, FUN = mean)
  pb$group_celltype <- interaction(pb$group, pb$celltype, drop = TRUE)
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(pb, ggplot2::aes(x = .data$group_celltype, y = .data$value, fill = .data$group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = alpha) +
    ggplot2::geom_jitter(width = 0.12, size = point_size, alpha = min(1, alpha + 0.15)) +
    scale_gleam_fill(palette = palette, continuous = FALSE) +
    ggplot2::labs(title = paste("Pseudobulk:", signature), x = "Group-Celltype", y = "Mean signature score") +
    do.call(gleam_theme, tp)
}
