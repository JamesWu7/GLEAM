#' Plot multiple signatures on spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param pathways Character vector of signature names (legacy argument name).
#' @param coords Spatial coordinates data.frame with x/y.
#' @param object Optional Seurat spatial object for in-slice plotting.
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
  coords = NULL,
  object = NULL,
  palette = "gleam_continuous",
  point_size = 1.2,
  alpha = 0.9,
  theme_params = list()
) {
  check_score_object(score)
  pathways <- intersect(pathways, rownames(score$score))
  if (length(pathways) < 1L) stop("No valid signatures provided.", call. = FALSE)
  tp <- resolve_text_params(theme_params)

  # Prefer native Seurat in-slice visualization for spatial Seurat objects.
  if (is.null(coords) && !is.null(object) && is_seurat_object(object) && requireNamespace("Seurat", quietly = TRUE)) {
    cells <- intersect(colnames(object), colnames(score$score))
    if (length(cells) < 1L) {
      stop("No overlapping cells between `object` and score matrix.", call. = FALSE)
    }

    feat_names <- character(length(pathways))
    for (i in seq_along(pathways)) {
      pw <- pathways[[i]]
      sig_col <- paste0("GLEAM_signature_", gsub("[^A-Za-z0-9_]+", "_", pw))
      vals <- rep(NA_real_, length(colnames(object)))
      names(vals) <- colnames(object)
      vals[cells] <- as.numeric(score$score[pw, cells])
      object[[sig_col]] <- vals
      feat_names[[i]] <- sig_col
    }

    p <- tryCatch({
      cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 128, continuous = TRUE) else palette
      img <- tryCatch(names(object@images)[[1]], error = function(e) NULL)
      Seurat::SpatialFeaturePlot(
        object = object,
        features = feat_names,
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

  ggplot2::ggplot(long, ggplot2::aes(x = .data$x, y = .data$y, color = .data$value)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    scale_gleam_color(palette = palette, continuous = TRUE) +
    ggplot2::facet_wrap(~ pathway) +
    ggplot2::coord_fixed() +
    ggplot2::scale_y_reverse() +
    ggplot2::labs(title = "Spatial multi-signature map", x = NULL, y = NULL, color = "Signature score") +
    do.call(gleam_theme, tp)
}
