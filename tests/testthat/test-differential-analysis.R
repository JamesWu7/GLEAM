test_that("pseudobulk level works", {
  data("pbmc_medium_matrix", package = "GLEAM")
  data("pbmc_medium_meta", package = "GLEAM")

  sc <- score_signature(
    expr = pbmc_medium_matrix,
    meta = pbmc_medium_meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  res <- test_signature(
    score = sc,
    group = "group",
    sample = "sample",
    celltype = "celltype",
    level = "pseudobulk",
    method = "wilcox",
    verbose = FALSE
  )

  expect_s3_class(res, "gleam_test")
  expect_true(nrow(res$table) > 0)
})
