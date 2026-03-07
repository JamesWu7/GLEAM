test_that("geneset supports data.frame and gmt path", {
  data("toy_expr", package = "GLEAM")

  gs_df <- data.frame(
    pathway = c("P_IFN", "P_IFN", "P_CYT", "P_CYT", "P_CYT"),
    gene = c("STAT1", "IRF1", "NKG7", "PRF1", "GZMB"),
    stringsAsFactors = FALSE
  )

  sc1 <- suppressWarnings(score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = gs_df,
    seurat = FALSE,
    method = "zmean",
    min_genes = 2,
    verbose = FALSE
  ))
  expect_s3_class(sc1, "gleam_score")

  gmt <- system.file("extdata", "genesets", "immune_small_example.gmt", package = "GLEAM")
  sc2 <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = gmt,
    seurat = FALSE,
    method = "auc",
    min_genes = 3,
    verbose = FALSE
  )
  expect_s3_class(sc2, "gleam_score")
})

test_that("ensemble scoring works", {
  data("toy_expr", package = "GLEAM")

  sc <- score_pathway(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "ensemble",
    min_genes = 3,
    verbose = FALSE
  )

  expect_s3_class(sc, "gleam_score")
  expect_true(all(dim(sc$score) > 0))
})
