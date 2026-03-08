test_that("aggregate_signature supports mean and fraction", {
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

  m1 <- aggregate_signature(sc, by = c("group", "celltype"), fun = "mean", long = FALSE)
  expect_true(is.matrix(m1))
  expect_true(nrow(m1) > 0)

  m2 <- aggregate_signature(sc, by = "group", fun = "fraction", threshold = 0.5, long = TRUE)
  expect_true(is.data.frame(m2))
  expect_true(all(c("pathway", "group_key", "value") %in% colnames(m2)))
})

test_that("test_signature supports sample and sample_celltype levels", {
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

  s1 <- test_signature(sc, group = "group", sample = "sample", level = "sample", verbose = FALSE)
  expect_s3_class(s1, "gleam_test")
  expect_true(nrow(s1$table) > 0)

  s2 <- test_signature(
    sc,
    group = "group",
    sample = "sample",
    celltype = "celltype",
    level = "sample_celltype",
    verbose = FALSE
  )
  expect_s3_class(s2, "gleam_test")
  expect_true(nrow(s2$table) > 0)
})
