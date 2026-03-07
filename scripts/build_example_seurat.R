if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
  stop("Please install Seurat/SeuratObject first.")
}

load("data/pbmc_medium_matrix.rda")
load("data/pbmc_medium_meta.rda")

if (requireNamespace("Seurat", quietly = TRUE)) {
  seu <- Seurat::CreateSeuratObject(counts = pbmc_medium_matrix, meta.data = pbmc_medium_meta)
  if ("sample" %in% colnames(pbmc_medium_meta)) {
    SeuratObject::DefaultAssay(seu) <- "RNA"
  }
  saveRDS(seu, file = "inst/extdata/pbmc_medium_seurat.rds")
  message("Saved inst/extdata/pbmc_medium_seurat.rds")
} else {
  stop("Seurat package is required for CreateSeuratObject.")
}
