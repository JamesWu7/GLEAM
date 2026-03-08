#' Plot signature score on spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param pathway Signature name (legacy argument name).
#' @param coords Spatial coordinates data.frame with x/y.
#' @param object Optional Seurat spatial object for in-slice plotting.
#' @param image Optional raster/image object. If provided, used as background.
#' @param palette Continuous palette name or colors.
#' @param split.by Optional metadata split variable.
#' @param point_size Point size.
#' @param alpha Point alpha.
#' @param theme_params Optional list passed to [gleam_theme()].
#'
#' @return A `ggplot` object.
#' @export
plot_spatial_score <- function(
  score,
  pathway,
  coords = NULL,
  object = NULL,
  image = NULL,
  palette = "gleam_continuous",
  split.by = NULL,
  point_size = 1.4,
  alpha = 0.9,
  theme_params = list()
) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`signature` not found in score matrix.", call. = FALSE)
  tp <- resolve_text_params(theme_params)

  # Prefer native Seurat spatial rendering when a spatial Seurat object is available.
  if (is.null(coords) && !is.null(object) && is_seurat_object(object) && requireNamespace("Seurat", quietly = TRUE)) {
    sig_col <- paste0("GLEAM_signature_", gsub("[^A-Za-z0-9_]+", "_", pathway))
    cells <- intersect(colnames(object), colnames(score$score))
    if (length(cells) < 1L) {
      stop("No overlapping cells between `object` and score matrix.", call. = FALSE)
    }
    vals <- rep(NA_real_, length(colnames(object)))
    names(vals) <- colnames(object)
    vals[cells] <- as.numeric(score$score[pathway, cells])
    object[[sig_col]] <- vals

    p <- tryCatch({
      cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 128, continuous = TRUE) else palette
      img <- tryCatch(names(object@images)[[1]], error = function(e) NULL)
      Seurat::SpatialFeaturePlot(
        object = object,
        features = sig_col,
        images = img,
        pt.size.factor = point_size,
        alpha = c(0.05, alpha),
        cols = cols
      )
    }, error = function(e) {
      warning(
        sprintf("Seurat SpatialFeaturePlot failed (%s). Falling back to coordinate scatter mode.", e$message),
        call. = FALSE
      )
      NULL
    })
    if (!is.null(p)) {
      if ("patchwork" %in% class(p) && requireNamespace("patchwork", quietly = TRUE)) {
        p <- p & do.call(gleam_theme, tp)
      } else if (inherits(p, "ggplot")) {
        p <- p + do.call(gleam_theme, tp)
      }
      return(p)
    }

    md <- tryCatch(extract_meta(object = object, seurat = TRUE), error = function(e) NULL)
    if (!is.null(md)) {
      if (all(c("x", "y") %in% colnames(md))) {
        coords <- data.frame(x = md$x, y = md$y, row.names = rownames(md))
      } else if (all(c("imagecol", "imagerow") %in% colnames(md))) {
        coords <- data.frame(x = md$imagecol, y = md$imagerow, row.names = rownames(md))
      } else if (all(c("col", "row") %in% colnames(md))) {
        coords <- data.frame(x = md$col, y = md$row, row.names = rownames(md))
      }
    }
  }

  if (is.null(coords)) {
    stop("Provide `coords`, or supply a Seurat spatial `object` for in-slice plotting.", call. = FALSE)
  }
  dat <- join_score_spatial(score, coords)
  dat$pathway_score <- as.numeric(score$score[pathway, dat$cell_id])

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$pathway_score))
  if (!is.null(image)) {
    p <- p + ggplot2::annotation_raster(image, xmin = min(dat$x), xmax = max(dat$x), ymin = min(dat$y), ymax = max(dat$y))
  }

  p <- p +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = paste("Spatial signature:", pathway), x = NULL, y = NULL, color = "Signature score") +
    ggplot2::scale_y_reverse() +
    do.call(gleam_theme, tp)

  if (!is.null(split.by)) {
    sp <- resolve_meta_var(score$meta, split.by, "split.by")
    dat$split <- as.factor(sp)
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data$x, y = .data$y, color = .data$pathway_score))
    if (!is.null(image)) {
      p <- p + ggplot2::annotation_raster(image, xmin = min(dat$x), xmax = max(dat$x), ymin = min(dat$y), ymax = max(dat$y))
    }
    p <- p +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      scale_gleam_color(palette = palette, continuous = TRUE) +
      ggplot2::coord_fixed() +
      ggplot2::facet_wrap(~ split) +
      ggplot2::labs(title = paste("Spatial signature:", pathway), x = NULL, y = NULL, color = "Signature score") +
      ggplot2::scale_y_reverse() +
      do.call(gleam_theme, tp)
  }

  p
}
