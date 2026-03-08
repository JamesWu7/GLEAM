test_that("Monocle3 dependency failure is clear when backend is requested", {
  if (requireNamespace("monocle3", quietly = TRUE)) {
    skip("monocle3 installed; missing-dependency behavior not applicable")
  }

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
    "Monocle3 is optional"
  )
})
