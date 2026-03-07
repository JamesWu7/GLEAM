#' Plot pathway score on spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param coords Spatial coordinates data.frame with x/y.
#' @param image Optional raster/image object. If provided, used as background.
#' @param palette Continuous palette name or colors.
#' @param split.by Optional metadata split variable.
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_score <- function(score, pathway, coords, image = NULL, palette = "gleam_continuous", split.by = NULL) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  dat <- join_score_spatial(score, coords)
  dat$pathway_score <- as.numeric(score$score[pathway, dat$cell_id])

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = pathway_score))
  if (!is.null(image)) {
    p <- p + ggplot2::annotation_raster(image, xmin = min(dat$x), xmax = max(dat$x), ymin = min(dat$y), ymax = max(dat$y))
  }

  p <- p +
    ggplot2::geom_point(size = 1.4, alpha = 0.9) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = paste("Spatial score:", pathway), x = "x", y = "y", color = pathway) +
    .theme_gleam()

  if (!is.null(split.by)) {
    sp <- resolve_meta_var(score$meta, split.by, "split.by")
    dat$split <- as.factor(sp)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, color = pathway_score)) +
      ggplot2::geom_point(size = 1.4, alpha = 0.9) +
      scale_gleam_color(palette = palette, continuous = TRUE) +
      ggplot2::coord_fixed() +
      ggplot2::facet_wrap(~ split) +
      ggplot2::labs(title = paste("Spatial score:", pathway), x = "x", y = "y", color = pathway) +
      .theme_gleam()
  }

  p
}
