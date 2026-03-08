#' Plot multiple pathways on spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param pathways Character vector of pathway names.
#' @param coords Spatial coordinates data.frame with x/y.
#' @param palette Continuous palette.
#' @param point_size Point size.
#' @param alpha Point alpha.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_multi <- function(
  score,
  pathways,
  coords,
  palette = "gleam_continuous",
  point_size = 1.2,
  alpha = 0.9,
  theme_params = list()
) {
  check_score_object(score)
  pathways <- intersect(pathways, rownames(score$score))
  if (length(pathways) < 1L) stop("No valid pathways provided.", call. = FALSE)

  dat <- join_score_spatial(score, coords)
  long <- do.call(rbind, lapply(pathways, function(pw) {
    data.frame(
      cell_id = dat$cell_id,
      x = dat$x,
      y = dat$y,
      pathway = pw,
      value = as.numeric(score$score[pw, dat$cell_id]),
      stringsAsFactors = FALSE
    )
  }))
  tp <- resolve_text_params(theme_params)

  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Spatial multi-pathway map", x = "x", y = "y") +
    do.call(gleam_theme, tp)
}
