#' Check Seurat compatibility mode
#'
#' @param object Seurat object.
#'
#' @return Character string describing detected Seurat mode.
#' @export
seurat_mode <- function(object) {
  if (!is_seurat_object(object)) return("not_seurat")
  v <- detect_seurat_version(object)
  if (is.na(v)) return("seurat_unknown")
  if (v >= 5L) return("seurat_v5")
  if (v == 4L) return("seurat_v4")
  sprintf("seurat_v%d", v)
}
