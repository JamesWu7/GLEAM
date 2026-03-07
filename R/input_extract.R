#' Resolve Seurat assay/layer with v5 compatibility
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

  # Default preference for v5: data, then counts.
  try_layers <- c("data", "counts")
  found <- NULL
  for (ly in try_layers) {
    ok <- tryCatch({
      SeuratObject::LayerData(object = object, assay = assay, layer = ly)
      TRUE
    }, error = function(e) FALSE)
    if (ok) {
      found <- ly
      break
    }
  }

  if (is.null(found) && !is.null(slot)) {
    return(list(assay = assay, layer = NULL, slot = slot))
  }
  if (is.null(found)) {
    found <- "data"
  }

  list(assay = assay, layer = found, slot = slot)
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

  has_seurat <- requireNamespace("SeuratObject", quietly = TRUE) || requireNamespace("Seurat", quietly = TRUE)
  if (!has_seurat) {
    stop("Seurat mode requested but Seurat/SeuratObject is not installed.", call. = FALSE)
  }

  info <- resolve_seurat_layer(object, assay = assay, layer = layer, slot = slot)

  mat <- tryCatch(
    {
      if (!is.null(info$layer)) {
        SeuratObject::LayerData(object = object, assay = info$assay, layer = info$layer)
      } else {
        SeuratObject::GetAssayData(object = object, assay = info$assay, slot = info$slot %||% "data")
      }
    },
    error = function(e) {
      if (!is.null(info$slot)) {
        SeuratObject::GetAssayData(object = object, assay = info$assay, slot = info$slot)
      } else {
        stop(sprintf("Failed to extract expression from Seurat object: %s", e$message), call. = FALSE)
      }
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
