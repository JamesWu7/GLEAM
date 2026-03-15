test_that("spatial score join and plotting work", {
  data("spatial_medium_expr", package = "GLEAM")
  data("spatial_medium_meta", package = "GLEAM")
  data("spatial_medium_coords", package = "GLEAM")

  sc <- score_signature(expr = spatial_medium_expr, meta = spatial_medium_meta, geneset = "immune_small", seurat = FALSE, method = "rank", min_genes = 3, verbose = FALSE)
  js <- join_score_spatial(sc, coords = spatial_medium_coords)
  expect_true(all(c("cell_id", "x", "y") %in% colnames(js)))

  p <- plot_spatial_score(sc, signature = rownames(sc$score)[1], coords = spatial_medium_coords)
  expect_s3_class(p, "ggplot")

  d <- test_signature_spatial(sc, region = "region", group = "group", sample = "sample", level = "sample_region")
  expect_s3_class(d, "gleam_test")
})
