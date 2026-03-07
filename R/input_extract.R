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

  # Seurat path: try image coordinates first, then metadata fallbacks.
  out <- tryCatch({
    imgs <- names(object@images)
    if (length(imgs) < 1) stop("no image")
    cdat <- SeuratObject::GetTissueCoordinates(object = object, image = imgs[[1]])
    cdat <- as.data.frame(cdat)
    if (!all(c("x", "y") %in% colnames(cdat))) {
      if (all(c("imagerow", "imagecol") %in% colnames(cdat))) {
        cdat <- data.frame(x = cdat$imagecol, y = cdat$imagerow, row.names = rownames(cdat))
      }
    } else {
      cdat <- cdat[, c("x", "y"), drop = FALSE]
    }
    cdat
  }, error = function(e) NULL)

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
