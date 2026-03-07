#' Extract embedding matrix
#'
#' @param object Seurat object.
#' @param reduction Reduction name.
#'
#' @return Embedding matrix.
#' @export
extract_embedding <- function(object, reduction = "umap") {
  extract_reduction(object = object, reduction = reduction)
}
