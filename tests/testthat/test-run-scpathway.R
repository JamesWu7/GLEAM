test_that("run_gleam end-to-end works", {
  data("toy_expr", package = "GLEAM")

  out <- run_gleam(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    group = "group",
    sample = "sample",
    celltype = "celltype",
    comparison = "within_celltype",
    target_celltype = "CD8_T",
    level = "sample",
    verbose = FALSE
  )

  expect_true(all(c("score", "test", "top_table") %in% names(out)))
  expect_s3_class(out$score, "gleam_score")
  expect_s3_class(out$test, "gleam_test")
})
