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
  sp_payload <- NULL
  if (is.null(coords) && !is.null(object) && is_seurat_object(object)) {
    sp_payload <- tryCatch(.extract_seurat_spatial_payload(object), error = function(e) NULL)
  }

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

    cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 128, continuous = TRUE) else palette
    img <- tryCatch({
      if (!is.null(sp_payload$image_name)) sp_payload$image_name else names(object@images)[[1]]
    }, error = function(e) NULL)
    sf_fun <- get0(".spatial_featureplot_compat", mode = "function", inherits = TRUE)
    if (is.null(sf_fun)) {
      sf <- list(plot = NULL, error = "internal helper .spatial_featureplot_compat() not found")
    } else {
      sf_args <- list(
        object = object,
        features = sig_col,
        image_name = img,
        point_size = point_size,
        alpha = alpha,
        cols = cols
      )
      fm <- names(formals(sf_fun))
      if (!("..." %in% fm)) {
        sf_args <- sf_args[intersect(names(sf_args), fm)]
      }
      sf <- tryCatch(
        do.call(sf_fun, sf_args),
        error = function(e) list(plot = NULL, error = conditionMessage(e))
      )
      if (!is.list(sf)) {
        sf <- list(plot = NULL, error = "invalid return from .spatial_featureplot_compat()")
      }
      if (is.null(sf$plot)) sf$plot <- NULL
      if (is.null(sf$error)) sf$error <- NULL
    }
    p <- sf$plot
    if (!is.null(p)) {
      if ("patchwork" %in% class(p) && requireNamespace("patchwork", quietly = TRUE)) {
        p <- p & do.call(gleam_theme, tp)
      } else if (inherits(p, "ggplot")) {
        p <- p + do.call(gleam_theme, tp)
      }
      return(p)
    }
    warning(
      sprintf("Seurat SpatialFeaturePlot failed (%s). Falling back to coordinate+tissue overlay mode.", sf$error %||% "unknown error"),
      call. = FALSE
    )
  }
  if (is.null(coords) && !is.null(sp_payload) && !is.null(sp_payload$coords)) {
    coords <- sp_payload$coords
  }
  if (is.null(coords) && is.data.frame(score$meta)) {
    coords <- tryCatch(.normalize_spatial_xy(score$meta), error = function(e) NULL)
  }
  if (is.null(image) && !is.null(sp_payload) && !is.null(sp_payload$image)) {
    image <- sp_payload$image
  }

  if (is.null(coords)) {
    stop("Provide `coords`, or supply a Seurat spatial `object` for in-slice plotting.", call. = FALSE)
  }
  dat <- join_score_spatial(score, coords)
  dat$pathway_score <- as.numeric(score$score[pathway, dat$cell_id])
  dat <- dat[order(dat$pathway_score), , drop = FALSE]

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
    do.call(gleam_theme, tp) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())

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
      do.call(gleam_theme, tp) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  p
}
