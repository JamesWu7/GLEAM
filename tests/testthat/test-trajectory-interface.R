test_that("trajectory mapping and testing work", {
  data("toy_expr", package = "GLEAM")
  sc <- score_signature(expr = toy_expr$expr, meta = toy_expr$meta, geneset = "immune_small", seurat = FALSE, method = "rank", min_genes = 3, verbose = FALSE)

  tr <- as_trajectory_data(sc, pseudotime = "pseudotime", lineage = "lineage")
  expect_true(all(c("cell_id", "pseudotime", "lineage") %in% colnames(tr)))

  tt <- test_signature_trajectory(sc, pseudotime = "pseudotime", lineage = "lineage", method = "spearman", verbose = FALSE)
  expect_s3_class(tt, "gleam_test")

  p <- plot_pseudotime_score(sc, pathway = rownames(sc$score)[1], pseudotime = "pseudotime", lineage = "lineage")
  expect_s3_class(p, "ggplot")
})
