test_that("plot_violin returns ggplot object", {
  skip_if_not_installed("ggplot2")

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

  p <- plot_violin(sc, signature = rownames(sc$score)[1], group = "group")
  expect_s3_class(p, "ggplot")
})
