#' @keywords internal
.spatial_featureplot_compat <- function(object, features, image_name = NULL, point_size = 1.4, alpha = 0.9, cols = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    return(list(plot = NULL, error = "Seurat not available"))
  }

  base_images <- list()
  if (!is.null(image_name) && is.character(image_name) && length(image_name) == 1L && nzchar(image_name)) {
    base_images <- list(
      list(images = image_name),
      list(image = image_name)
    )
  } else {
    base_images <- list(list(), list(images = NULL), list(image = NULL))
  }

  style_args <- list(
    list(pt.size.factor = point_size, alpha = c(0.05, alpha), cols = cols),
    list(pt.size.factor = point_size, alpha = c(0.05, alpha)),
    list(pt.size.factor = point_size, cols = cols),
    list(alpha = c(0.05, alpha), cols = cols),
    list(pt.size.factor = point_size),
    list(alpha = c(0.05, alpha)),
    list(cols = cols),
    list()
  )

  last_err <- NULL
  for (img_arg in base_images) {
    for (sty in style_args) {
      args <- c(list(object = object, features = features), img_arg, sty)
      args <- args[!vapply(args, is.null, logical(1))]
      out <- tryCatch(
        do.call(Seurat::SpatialFeaturePlot, args),
        error = function(e) e
      )
      if (!inherits(out, "error")) {
        return(list(plot = out, error = NULL))
      }
      last_err <- conditionMessage(out)
    }
  }
  list(plot = NULL, error = last_err %||% "SpatialFeaturePlot failed")
}
