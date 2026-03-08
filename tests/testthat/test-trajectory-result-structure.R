test_that("trajectory result structure is attached and stable", {
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

  tr <- GLEAM:::as_trajectory_result(
    score = sc,
    pseudotime = "pseudotime",
    lineage = "lineage",
    backend = "internal"
  )
  expect_s3_class(tr, "gleam_trajectory_result")
  expect_true(all(c("backend", "pseudotime", "lineage", "meta", "params", "version_info") %in% names(tr)))
  expect_equal(length(tr$pseudotime), ncol(sc$score))

  tt <- test_signature_trajectory(
    score = sc,
    pseudotime = "pseudotime",
    lineage = "lineage",
    method = "spearman",
    backend = "internal",
    verbose = FALSE
  )
  expect_s3_class(tt, "gleam_test")
  expect_true("trajectory" %in% names(tt$comparison))
  expect_equal(tt$comparison$trajectory$backend, "internal")
})
