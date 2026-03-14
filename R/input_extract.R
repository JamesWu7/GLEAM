#' Resolve Seurat assay/layer with v4/v5 compatibility
#'
#' @param object Seurat object.
#' @param assay Assay name.
#' @param layer Layer name.
#' @param slot Legacy slot name.
#'
#' @return List with assay and layer/slot.
#' @keywords internal
resolve_seurat_layer <- function(object, assay = NULL, layer = NULL, slot = NULL) {
  assay <- assay %||% tryCatch(SeuratObject::DefaultAssay(object), error = function(e) NULL)
  if (is.null(assay)) {
    stop("Unable to resolve assay from Seurat object; please supply `assay`.", call. = FALSE)
  }

  if (!is.null(layer)) {
    return(list(assay = assay, layer = layer, slot = slot))
  }

  version_major <- detect_seurat_version(object)

  # v5 prefers layers; v4 primarily uses slots.
  if (!is.na(version_major) && version_major >= 5L) {
    try_layers <- c("data", "counts")
    for (ly in try_layers) {
      ok <- tryCatch({
        SeuratObject::LayerData(object = object, assay = assay, layer = ly)
        TRUE
      }, error = function(e) FALSE)
      if (ok) {
        return(list(assay = assay, layer = ly, slot = slot))
      }
    }
  }

  list(assay = assay, layer = NULL, slot = slot %||% "data")
}

#' @keywords internal
.is_valid_expr_matrix <- function(x) {
  if (is.null(x)) return(FALSE)
  ok_cls <- is.matrix(x) || inherits(x, "dgCMatrix")
  if (!ok_cls) return(FALSE)
  if (nrow(x) < 1L || ncol(x) < 1L) return(FALSE)
  !is.null(rownames(x)) && !is.null(colnames(x))
}

#' Extract expression matrix from input
#'
#' @param object Seurat object.
#' @param expr Matrix input when `seurat = FALSE`.
#' @param assay Assay name.
#' @param layer Layer name.
#' @param slot Legacy slot.
#' @param seurat Input mode.
#'
#' @return Matrix or `dgCMatrix`.
#' @keywords internal
extract_expr <- function(object = NULL, expr = NULL, assay = NULL, layer = NULL, slot = NULL, seurat = TRUE) {
  if (!seurat) {
    return(as_expr_matrix(expr))
  }

  info <- resolve_seurat_layer(object, assay = assay, layer = layer, slot = slot)

  mat <- tryCatch(
    {
      if (!is.null(info$layer)) {
        SeuratObject::LayerData(object = object, assay = info$assay, layer = info$layer)
      } else {
        SeuratObject::GetAssayData(object = object, assay = info$assay, slot = info$slot)
      }
    },
    error = function(e) {
      # Fallback path for mixed Seurat v4/v5 environments.
      tryCatch(
        SeuratObject::GetAssayData(object = object, assay = info$assay, slot = info$slot %||% "data"),
        error = function(e2) {
          stop(sprintf("Failed to extract expression from Seurat object: %s", e2$message), call. = FALSE)
        }
      )
    }
  )

  if (!.is_valid_expr_matrix(mat)) {
    fallback_expr <- NULL
    alt_layers <- unique(c("data", "counts"))
    if (!is.null(info$layer)) {
      alt_layers <- c(info$layer, setdiff(alt_layers, info$layer))
    }
    for (ly in alt_layers) {
      cand <- tryCatch(
        SeuratObject::LayerData(object = object, assay = info$assay, layer = ly),
        error = function(e) NULL
      )
      if (.is_valid_expr_matrix(cand)) {
        fallback_expr <- cand
        break
      }
    }
    if (is.null(fallback_expr)) {
      for (sl in c("data", "counts")) {
        cand <- tryCatch(
          SeuratObject::GetAssayData(object = object, assay = info$assay, slot = sl),
          error = function(e) NULL
        )
        if (.is_valid_expr_matrix(cand)) {
          fallback_expr <- cand
          break
        }
      }
    }
    if (!is.null(fallback_expr)) {
      mat <- fallback_expr
    }
  }

  if (!.is_valid_expr_matrix(mat)) {
    stop(
      "Failed to extract a non-empty expression matrix with row/column names from Seurat object.",
      " Try specifying `assay`, `layer = \"counts\"`, or `slot = \"counts\"`.",
      call. = FALSE
    )
  }

  as_expr_matrix(mat)
}

#' Extract metadata from input
#'
#' @param object Seurat object.
#' @param meta Metadata data.frame for matrix mode.
#' @param expr Expression matrix for matrix mode.
#' @param seurat Input mode.
#'
#' @return Metadata data.frame aligned to cells.
#' @keywords internal
extract_meta <- function(object = NULL, meta = NULL, expr = NULL, seurat = TRUE) {
  if (!seurat) {
    if (is.null(meta)) {
      meta <- data.frame(cell_id = colnames(expr), row.names = colnames(expr), stringsAsFactors = FALSE)
    }
    if (is.null(rownames(meta))) {
      rownames(meta) <- colnames(expr)
    }
    return(meta)
  }

  md <- tryCatch(
    as.data.frame(object[[]], stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(md)) {
    stop("Failed to extract metadata from Seurat object.", call. = FALSE)
  }
  if (is.null(rownames(md))) {
    rownames(md) <- colnames(object)
  }
  md
}

#' Extract reduction embedding matrix
#'
#' @param object Seurat object.
#' @param reduction Reduction name (for example: `umap`, `pca`).
#'
#' @return Numeric embedding matrix.
#' @keywords internal
extract_reduction <- function(object, reduction = "umap") {
  if (!is_seurat_object(object)) {
    stop("`object` must be Seurat for reduction extraction.", call. = FALSE)
  }

  emb <- tryCatch(
    SeuratObject::Embeddings(object = object, reduction = reduction),
    error = function(e) NULL
  )
  if (is.null(emb)) {
    stop(sprintf("Reduction '%s' not available in Seurat object.", reduction), call. = FALSE)
  }
  emb
}

#' Extract embeddings from generic input
#'
#' @param object Seurat object.
#' @param embeddings Numeric matrix, optional direct input.
#' @param reduction Reduction name when extracting from Seurat.
#' @param seurat Input mode.
#'
#' @return Embedding matrix.
#' @keywords internal
extract_embeddings <- function(object = NULL, embeddings = NULL, reduction = "umap", seurat = TRUE) {
  if (!is.null(embeddings)) {
    return(as.matrix(embeddings))
  }
  if (!seurat) {
    stop("Provide `embeddings` matrix when `seurat = FALSE`.", call. = FALSE)
  }
  extract_reduction(object = object, reduction = reduction)
}

#' Check whether input contains spatial information
#'
#' @param object Seurat object or list.
#' @param meta Metadata data.frame.
#'
#' @return Logical.
#' @keywords internal
is_spatial_object <- function(object = NULL, meta = NULL) {
  if (is_seurat_object(object)) {
    has_images <- tryCatch(length(object@images) > 0, error = function(e) FALSE)
    if (has_images) return(TRUE)
  }
  if (!is.null(meta) && is.data.frame(meta)) {
    xy <- c("x", "y")
    rc <- c("row", "col")
    return(all(xy %in% colnames(meta)) || all(rc %in% colnames(meta)))
  }
  FALSE
}

#' Extract spatial coordinates
#'
#' @param object Seurat object.
#' @param meta Metadata for matrix mode.
#' @param coords Optional direct coordinate data.frame.
#' @param seurat Input mode.
#'
#' @return Data.frame with columns `x` and `y`.
#' @keywords internal
extract_spatial_coords <- function(object = NULL, meta = NULL, coords = NULL, seurat = TRUE) {
  if (!is.null(coords)) {
    coords <- as.data.frame(coords)
    check_required_columns(coords, c("x", "y"))
    return(coords[, c("x", "y"), drop = FALSE])
  }

  if (!seurat) {
    if (is.null(meta)) stop("`meta` required to infer spatial coordinates.", call. = FALSE)
    meta <- as.data.frame(meta)
    if (all(c("x", "y") %in% colnames(meta))) {
      return(meta[, c("x", "y"), drop = FALSE])
    }
    if (all(c("row", "col") %in% colnames(meta))) {
      out <- meta[, c("col", "row"), drop = FALSE]
      colnames(out) <- c("x", "y")
      return(out)
    }
    stop("No spatial coordinates found. Provide `coords` or `meta` columns x/y or row/col.", call. = FALSE)
  }

  # Seurat path: try unified spatial payload first, then metadata fallbacks.
  payload <- tryCatch(.extract_seurat_spatial_payload(object), error = function(e) NULL)
  out <- if (!is.null(payload)) payload$coords else NULL

  if (!is.null(out)) return(out)

  md <- extract_meta(object = object, seurat = TRUE)
  extract_spatial_coords(meta = md, seurat = FALSE)
}

#' Extract spatial metadata
#'
#' @param object Seurat object.
#' @param meta Metadata for matrix mode.
#' @param seurat Input mode.
#'
#' @return Metadata data.frame.
#' @keywords internal
extract_spatial_meta <- function(object = NULL, meta = NULL, seurat = TRUE) {
  if (seurat) {
    return(extract_meta(object = object, seurat = TRUE))
  }
  if (is.null(meta)) {
    stop("`meta` must be provided for matrix-based spatial workflows.", call. = FALSE)
  }
  as.data.frame(meta, stringsAsFactors = FALSE)
}

#' @keywords internal
.normalize_spatial_xy <- function(df) {
  if (is.null(df)) return(NULL)
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (all(c("x", "y") %in% colnames(df))) {
    out <- df[, c("x", "y"), drop = FALSE]
  } else if (all(c("coords_x", "coords_y") %in% colnames(df))) {
    out <- data.frame(x = df$coords_x, y = df$coords_y, row.names = rownames(df))
  } else if (all(c("x_coord", "y_coord") %in% colnames(df))) {
    out <- data.frame(x = df$x_coord, y = df$y_coord, row.names = rownames(df))
  } else if (all(c("imagecol", "imagerow") %in% colnames(df))) {
    out <- data.frame(x = df$imagecol, y = df$imagerow, row.names = rownames(df))
  } else if (all(c("array_col", "array_row") %in% colnames(df))) {
    out <- data.frame(x = df$array_col, y = df$array_row, row.names = rownames(df))
  } else if (all(c("col", "row") %in% colnames(df))) {
    out <- data.frame(x = df$col, y = df$row, row.names = rownames(df))
  } else if (ncol(df) >= 2L) {
    out <- data.frame(x = df[[1]], y = df[[2]], row.names = rownames(df))
  } else {
    return(NULL)
  }
  out$x <- as.numeric(out$x)
  out$y <- as.numeric(out$y)
  out
}

#' @keywords internal
.as_raster_spatial_image <- function(img) {
  if (is.null(img)) return(NULL)
  if (inherits(img, "raster")) return(img)
  if (is.matrix(img)) return(grDevices::as.raster(img))

  if (is.array(img) && length(dim(img)) == 3L && dim(img)[3] >= 3L) {
    arr <- img[, , seq_len(3), drop = FALSE]
    arr <- suppressWarnings(as.numeric(arr))
    arr <- array(arr, dim = dim(img)[1:3])
    if (!all(is.na(arr)) && max(arr, na.rm = TRUE) > 1) {
      arr <- arr / 255
    }
    arr <- pmax(0, pmin(1, arr))
    rgb_mat <- grDevices::rgb(arr[, , 1], arr[, , 2], arr[, , 3])
    dim(rgb_mat) <- dim(arr)[1:2]
    return(grDevices::as.raster(rgb_mat))
  }
  NULL
}

#' @keywords internal
.extract_seurat_spatial_payload <- function(object, image = NULL) {
  if (!is_seurat_object(object)) return(NULL)
  img_names <- tryCatch(names(object@images), error = function(e) character())
  if (length(img_names) < 1L) return(NULL)

  image_name <- if (!is.null(image) && is.character(image) && length(image) == 1L && image %in% img_names) {
    image
  } else {
    img_names[[1]]
  }

  img_obj <- tryCatch(object@images[[image_name]], error = function(e) NULL)

  coords <- tryCatch({
    tc <- SeuratObject::GetTissueCoordinates(object = object, image = image_name)
    .normalize_spatial_xy(tc)
  }, error = function(e) NULL)
  if (is.null(coords) && requireNamespace("Seurat", quietly = TRUE)) {
    coords <- tryCatch({
      tc <- Seurat::GetTissueCoordinates(object = object, image = image_name)
      .normalize_spatial_xy(tc)
    }, error = function(e) NULL)
  }

  if (is.null(coords) && !is.null(img_obj)) {
    bounds <- tryCatch(attr(img_obj, "boundaries"), error = function(e) NULL)
    if (!is.null(bounds) && length(bounds) > 0L) {
      for (b in bounds) {
        cdat <- tryCatch(attr(b, "coords"), error = function(e) NULL)
        cells <- tryCatch(attr(b, "cells"), error = function(e) NULL)
        cdat <- .normalize_spatial_xy(cdat)
        if (!is.null(cdat)) {
          if (!is.null(cells) && length(cells) == nrow(cdat)) {
            rownames(cdat) <- as.character(cells)
          }
          coords <- cdat
          break
        }
      }
    }
  }

  if (is.null(coords)) {
    md <- tryCatch(extract_meta(object = object, seurat = TRUE), error = function(e) NULL)
    coords <- .normalize_spatial_xy(md)
  }

  bg <- NULL
  if (!is.null(img_obj)) {
    # Seurat spatial image storage differs across versions/classes:
    # attribute("image"), @image slot, or methods::slot(image, "image").
    bg <- tryCatch(attr(img_obj, "image"), error = function(e) NULL)
    if (is.null(bg)) {
      bg <- tryCatch(img_obj@image, error = function(e) NULL)
    }
    if (is.null(bg) && methods::isS4(img_obj) && "image" %in% methods::slotNames(img_obj)) {
      bg <- tryCatch(methods::slot(img_obj, "image"), error = function(e) NULL)
    }
    bg <- .as_raster_spatial_image(bg)
  }

  if (is.null(coords) && is.null(bg)) {
    return(NULL)
  }
  list(coords = coords, image = bg, image_name = image_name)
}
