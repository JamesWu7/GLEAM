test_that("score_signature stores canonical method names and provenance fields", {
  data("toy_expr", package = "GLEAM")

  sc <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "zmean",
    min_genes = 3,
    verbose = FALSE
  )

  expect_identical(sc$method, "zscore")
  expect_true("geneset_size_original" %in% names(sc$params))
  expect_true("geneset_size_matched" %in% names(sc$params))
  expect_true("genes_dropped" %in% names(sc$params))
  expect_true("expression_layer_used" %in% names(sc$params))
  expect_true("method_parameters_used" %in% names(sc$params))
})
