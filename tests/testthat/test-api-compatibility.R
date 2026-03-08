test_that("signature aliases are available and compatible with pathway APIs", {
  data("toy_expr", package = "GLEAM")

  sc_sig <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )
  sc_old <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  expect_s3_class(sc_sig, "gleam_score")
  expect_equal(dim(sc_sig$score), dim(sc_old$score))

  tt_sig <- test_signature(sc_sig, group = "group", level = "cell", verbose = FALSE)
  tt_old <- test_pathway(sc_old, group = "group", level = "cell", verbose = FALSE)
  expect_s3_class(tt_sig, "gleam_test")
  expect_equal(colnames(tt_sig$table), colnames(tt_old$table))
})

test_that("legacy wrappers remain available", {
  expect_true(is.function(compare_methods))
  expect_true(is.function(run_scpathway))
  expect_true(is.function(score_pathway))
  expect_true(is.function(test_pathway))
})
