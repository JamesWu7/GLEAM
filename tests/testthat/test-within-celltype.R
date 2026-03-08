test_that("within-celltype comparison runs", {
  data("toy_expr", package = "GLEAM")
  sc <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  cmp <- test_signature(
    score = sc,
    group = "group",
    sample = "sample",
    celltype = "celltype",
    target_celltype = "CD8_T",
    level = "sample",
    method = "wilcox",
    verbose = FALSE
  )

  expect_s3_class(cmp, "gleam_test")
  expect_true(nrow(cmp$table) > 0)
})
