test_that("packaged Seurat spatial objects expose valid coordinates", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  full_example_path <- function(file_name) {
    p <- system.file("extdata", "full_examples", file_name, package = "GLEAM")
    if (nzchar(p) && file.exists(p)) return(p)
    p <- file.path("inst", "extdata", "full_examples", file_name)
    if (file.exists(p)) return(p)
    NA_character_
  }

  files <- c("stxBrain_anterior1_seurat.rds", "stxBrain_posterior1_seurat.rds")
  for (f in files) {
    p <- full_example_path(f)
    skip_if(is.na(p), paste("missing spatial example:", f))

    obj <- readRDS(p)
    expect_true(GLEAM:::is_spatial_object(object = obj))

    coords <- GLEAM:::extract_spatial_coords(object = obj)
    expect_s3_class(coords, "data.frame")
    expect_true(all(c("x", "y") %in% colnames(coords)))
    expect_true(nrow(coords) >= ncol(obj))
    expect_true(all(colnames(obj) %in% rownames(coords)))
  }
})

test_that("coordinate alignment remains valid through scoring + spatial join", {
  skip_on_cran()
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")

  p <- system.file("extdata", "full_examples", "stxBrain_anterior1_seurat.rds", package = "GLEAM")
  if (!(nzchar(p) && file.exists(p))) {
    p <- file.path("inst", "extdata", "full_examples", "stxBrain_anterior1_seurat.rds")
  }
  skip_if_not(file.exists(p), "missing spatial Seurat example")

  obj <- readRDS(p)
  assay <- SeuratObject::DefaultAssay(obj)
  expr <- SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  cells <- colnames(expr)[seq_len(min(300L, ncol(expr)))]
  genes <- rownames(expr)[seq_len(min(60L, nrow(expr)))]
  expr <- expr[genes, cells, drop = FALSE]

  meta <- as.data.frame(obj[[]], stringsAsFactors = FALSE)
  meta <- meta[cells, , drop = FALSE]

  gs <- list(Spatial_signature_test = genes[seq_len(min(20L, length(genes)))])
  sc <- score_signature(
    expr = expr,
    meta = meta,
    geneset = gs,
    geneset_source = "list",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  coords <- GLEAM:::extract_spatial_coords(object = obj)
  coords <- coords[colnames(sc$score), c("x", "y"), drop = FALSE]
  js <- join_score_spatial(sc, coords = coords)

  expect_equal(nrow(js), ncol(sc$score))
  expect_true(all(js$cell_id == colnames(sc$score)))
  expect_false(any(is.na(js$x)))
  expect_false(any(is.na(js$y)))
})

test_that("extract_spatial_coords reports actionable diagnostics for non-spatial Seurat objects", {
  skip_on_cran()
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Matrix")

  mat <- Matrix::Matrix(matrix(
    rpois(120, lambda = 2),
    nrow = 12,
    ncol = 10,
    dimnames = list(paste0("g", seq_len(12)), paste0("c", seq_len(10)))
  ), sparse = TRUE)
  obj <- SeuratObject::CreateSeuratObject(counts = mat)

  expect_error(
    GLEAM:::extract_spatial_coords(object = obj),
    regexp = "No spatial images/FOV found|Provide explicit `coords`|Metadata fallback failed"
  )
})
