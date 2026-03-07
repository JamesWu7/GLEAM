#' Plot multiple pathways on spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param pathways Character vector of pathway names.
#' @param coords Spatial coordinates data.frame with x/y.
#' @param palette Continuous palette.
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_multi <- function(score, pathways, coords, palette = "gleam_continuous") {
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

  ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, color = value)) +
    ggplot2::geom_point(size = 1.2, alpha = 0.9) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = "Spatial multi-pathway map", x = "x", y = "y") +
    .theme_gleam()
}
