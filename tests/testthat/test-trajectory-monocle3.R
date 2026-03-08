test_that("Monocle3 backend path is exercised when installed", {
  skip_if_not_installed("monocle3")

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

  fake_cds <- structure(list(), class = "cell_data_set")
  expect_error(
    test_signature_trajectory(
      score = sc,
      pseudotime = fake_cds,
      lineage = fake_cds,
      backend = "monocle3",
      verbose = FALSE
    ),
    "Failed to extract pseudotime|cell_data_set"
  )
})
