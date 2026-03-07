test_that("rank scoring returns pathway by cell matrix", {
  data("toy_expr", package = "GLEAM")
  data("immune_small", package = "GLEAM")

  sc <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = immune_small,
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  expect_s3_class(sc, "gleam_score")
  expect_true(nrow(sc$score) > 0)
  expect_equal(ncol(sc$score), ncol(toy_expr$expr))
})
