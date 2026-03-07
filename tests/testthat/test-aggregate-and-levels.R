test_that("aggregate_pathway supports mean and fraction", {
  data("toy_expr", package = "scPathway")
  sc <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  m1 <- aggregate_pathway(sc, by = c("group", "celltype"), fun = "mean", long = FALSE)
  expect_true(is.matrix(m1))
  expect_true(nrow(m1) > 0)

  m2 <- aggregate_pathway(sc, by = "group", fun = "fraction", threshold = 0.5, long = TRUE)
  expect_true(is.data.frame(m2))
  expect_true(all(c("pathway", "group_key", "value") %in% colnames(m2)))
})

test_that("test_pathway supports sample and sample_celltype levels", {
  data("toy_expr", package = "scPathway")
  sc <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  s1 <- test_pathway(sc, group = "group", sample = "sample", level = "sample", verbose = FALSE)
  expect_s3_class(s1, "scpathway_test")
  expect_true(nrow(s1$table) > 0)

  s2 <- test_pathway(
    sc,
    group = "group",
    sample = "sample",
    celltype = "celltype",
    level = "sample_celltype",
    verbose = FALSE
  )
  expect_s3_class(s2, "scpathway_test")
  expect_true(nrow(s2$table) > 0)
})
