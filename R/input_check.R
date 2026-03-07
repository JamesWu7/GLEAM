#' Check whether object is Seurat-like
#'
#' @param object Input object.
#'
#' @return Logical.
#' @keywords internal
is_seurat_object <- function(object) {
  inherits(object, "Seurat") || (!is.null(object) && "Seurat" %in% class(object))
}

#' Detect Seurat major version
#'
#' @param object Seurat object.
#'
#' @return Integer major version (4, 5) when detectable, otherwise `NA_integer_`.
#' @keywords internal
detect_seurat_version <- function(object = NULL) {
  if (!is_seurat_object(object)) {
    return(NA_integer_)
  }

  v <- NA_character_
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    v <- as.character(utils::packageVersion("SeuratObject"))
  } else if (requireNamespace("Seurat", quietly = TRUE)) {
    v <- as.character(utils::packageVersion("Seurat"))
  }

  if (is.na(v) || !nzchar(v)) return(NA_integer_)
  as.integer(strsplit(v, "\\.")[[1]][1])
}

#' Validate input mode
#'
#' @param object Seurat object.
#' @param expr Expression matrix.
#' @param seurat Logical, whether Seurat mode is used.
#'
#' @return Invisibly TRUE.
#' @keywords internal
check_input_mode <- function(object, expr, seurat = TRUE) {
  if (!is.logical(seurat) || length(seurat) != 1L) {
    stop("`seurat` must be TRUE or FALSE.", call. = FALSE)
  }

  if (seurat) {
    if (is.null(object)) {
      stop("`object` must be provided when `seurat = TRUE`.", call. = FALSE)
    }
    if (!is_seurat_object(object)) {
      stop("`object` must be a Seurat object when `seurat = TRUE`.", call. = FALSE)
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat mode requested but Seurat/SeuratObject is not installed.", call. = FALSE)
    }
  } else {
    if (is.null(expr)) {
      stop("`expr` must be provided when `seurat = FALSE`.", call. = FALSE)
    }
    as_expr_matrix(expr)
  }

  invisible(TRUE)
}

#' Validate expression and metadata alignment
#'
#' @param expr Expression matrix.
#' @param meta Metadata data.frame.
#'
#' @return Invisibly TRUE.
#' @keywords internal
check_expr_meta <- function(expr, meta = NULL) {
  expr <- as_expr_matrix(expr)

  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("`expr` must have both rownames (genes) and colnames (cells).", call. = FALSE)
  }

  if (!is.null(meta)) {
    if (!is.data.frame(meta)) {
      stop("`meta` must be a data.frame.", call. = FALSE)
    }
    if (nrow(meta) != ncol(expr)) {
      stop("nrow(meta) must equal ncol(expr).", call. = FALSE)
    }
    if (is.null(rownames(meta))) {
      rownames(meta) <- colnames(expr)
    }
  }

  invisible(TRUE)
}
