# Derive lightweight Seurat subsets from full example objects.
# This script is for maintainers and is not run during package load/check.

full_dir <- file.path("inst", "extdata", "full_examples")
out_dir <- file.path("inst", "extdata", "test_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat/SeuratObject is required to derive subset objects.")
}

subset_n <- 500L

subset_and_save <- function(infile, outfile, n = subset_n) {
  obj <- readRDS(infile)
  if (!inherits(obj, "Seurat")) {
    stop(sprintf("Input is not a Seurat object: %s", infile))
  }
  cells <- colnames(obj)
  keep <- utils::head(cells, min(length(cells), n))
  sub <- obj[, keep]
  saveRDS(sub, outfile)
  message(sprintf("[GLEAM] wrote %s (%d cells)", outfile, ncol(sub)))
}

subset_and_save(
  infile = file.path(full_dir, "ifnb_seurat.rds"),
  outfile = file.path(out_dir, "ifnb_seurat_subset.rds")
)
subset_and_save(
  infile = file.path(full_dir, "stxBrain_anterior1_seurat.rds"),
  outfile = file.path(out_dir, "stxBrain_anterior1_seurat_subset.rds")
)
subset_and_save(
  infile = file.path(full_dir, "stxBrain_posterior1_seurat.rds"),
  outfile = file.path(out_dir, "stxBrain_posterior1_seurat_subset.rds")
)
