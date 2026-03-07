test_that("cell-level testing returns required fields", {
  data("toy_expr", package = "GLEAM")
  sc <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  tt <- test_pathway(sc, group = "group", level = "cell", method = "wilcox", verbose = FALSE)
  expect_s3_class(tt, "gleam_test")

  req <- c(
    "pathway", "comparison_type", "group1", "group2", "celltype", "level",
    "effect_size", "median_group1", "median_group2", "diff_median",
    "p_value", "p_adj", "n_group1", "n_group2"
  )
  expect_true(all(req %in% colnames(tt$table)))
})
