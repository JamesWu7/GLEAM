test_that("scoring method registry is available", {
  m <- list_scoring_methods()
  expect_true(is.data.frame(m))
  expect_true(all(c("method_name", "category", "dependency") %in% colnames(m)))
  expect_true("rank" %in% m$method_name)
  expect_true("ensemble" %in% m$method_name)
})

test_that("new native scoring methods run", {
  data("toy_expr", package = "GLEAM")
  meths <- c("mean", "scaled_mean", "robust", "singscore_like", "ssgsea_like")
  for (m in meths) {
    sc <- score_pathway(
      expr = toy_expr$expr,
      meta = toy_expr$meta,
      geneset = "immune_small",
      seurat = FALSE,
      method = m,
      min_genes = 3,
      verbose = FALSE
    )
    expect_s3_class(sc, "gleam_score")
  }
})
