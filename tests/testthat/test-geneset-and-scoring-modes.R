test_that("geneset supports data.frame and gmt path", {
  data("toy_expr", package = "GLEAM")

  gs_df <- data.frame(
    pathway = c("P_IFN", "P_IFN", "P_CYT", "P_CYT", "P_CYT"),
    gene = c("STAT1", "IRF1", "NKG7", "PRF1", "GZMB"),
    stringsAsFactors = FALSE
  )

  sc1 <- suppressWarnings(score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = gs_df,
    seurat = FALSE,
    method = "zscore",
    min_genes = 2,
    verbose = FALSE
  ))
  expect_s3_class(sc1, "gleam_score")

  gmt <- system.file("extdata", "genesets", "immune_small_example.gmt", package = "GLEAM")
  sc2 <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = gmt,
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )
  expect_s3_class(sc2, "gleam_score")
})

test_that("ensemble scoring works", {
  data("toy_expr", package = "GLEAM")

  sc <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "ensemble",
    ensemble_methods = c("rank", "mean", "zscore"),
    ensemble_standardize = "zscore",
    ensemble_combine = "mean",
    min_genes = 3,
    verbose = FALSE
  )

  expect_s3_class(sc, "gleam_score")
  expect_true(all(dim(sc$score) > 0))
  expect_true(is.numeric(sc$score))
  expect_true(is.list(sc$params$method_parameters_used))
  expect_identical(sc$params$method_parameters_used$ensemble_standardize, "zscore")
  expect_true(all(abs(rowMeans(sc$score, na.rm = TRUE)) < 1))

  sc_w <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "ensemble",
    ensemble_methods = c("rank", "mean", "zscore"),
    ensemble_standardize = "zscore",
    ensemble_weights = c(rank = 1.2, mean = 0.8, zscore = 1),
    min_genes = 3,
    verbose = FALSE
  )
  expect_s3_class(sc_w, "gleam_score")
  expect_true(is.numeric(sc_w$score))
})

test_that("signed up/down signatures are supported", {
  data("toy_expr", package = "GLEAM")
  sig <- list(
    IFN_signed = list(
      up = c("STAT1", "IRF1", "ISG15"),
      down = c("LYZ", "S100A8")
    )
  )

  sc <- score_signature(
    expr = toy_expr$expr,
    meta = toy_expr$meta,
    geneset = sig,
    geneset_source = "list",
    seurat = FALSE,
    method = "rank",
    min_genes = 2,
    verbose = FALSE
  )

  expect_s3_class(sc, "gleam_score")
  expect_equal(rownames(sc$score), "IFN_signed")
  expect_true(is.numeric(sc$score[1, ]))
  expect_true("geneset_size_original" %in% names(sc$params))
  expect_true("geneset_size_matched" %in% names(sc$params))
  expect_true("genes_dropped" %in% names(sc$params))
  expect_true("expression_layer_used" %in% names(sc$params))
  expect_true("method_parameters_used" %in% names(sc$params))
})

test_that("as_geneset preserves signed up/down structure", {
  sig <- list(
    IFN_signed = list(
      up = c("STAT1", "IRF1", "ISG15"),
      down = c("LYZ", "S100A8")
    )
  )
  gs <- as_geneset(sig)
  expect_true(is.list(gs$IFN_signed))
  expect_setequal(gs$IFN_signed$up, c("STAT1", "IRF1", "ISG15"))
  expect_setequal(gs$IFN_signed$down, c("LYZ", "S100A8"))
})
