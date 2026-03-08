files <- c(
  file.path("inst", "extdata", "full_examples", "ifnb_seurat.rds"),
  file.path("inst", "extdata", "full_examples", "stxBrain_anterior1_seurat.rds"),
  file.path("inst", "extdata", "full_examples", "stxBrain_posterior1_seurat.rds")
)

if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
  message("[GLEAM] Seurat/SeuratObject not installed; skipping full-object structural audit.")
  quit(save = "no", status = 0)
}

for (f in files) {
  if (!file.exists(f)) {
    warning("[GLEAM] Missing example object: ", f, call. = FALSE)
    next
  }
  obj <- readRDS(f)
  if (!inherits(obj, "Seurat")) {
    warning("[GLEAM] Not a Seurat object: ", f, call. = FALSE)
    next
  }

  reductions <- names(obj@reductions)
  md <- colnames(obj@meta.data)
  has_xy <- all(c("x", "y") %in% md)
  has_images <- tryCatch(length(obj@images) > 0, error = function(e) FALSE)

  message("[GLEAM] ", basename(f))
  message("  - dim: ", nrow(obj), " genes x ", ncol(obj), " cells")
  message("  - reductions: ", if (length(reductions) == 0) "<none>" else paste(reductions, collapse = ", "))
  message("  - has x/y metadata: ", has_xy)
  message("  - has image slots: ", has_images)
}
