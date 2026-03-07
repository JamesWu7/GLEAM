test_that("export and method comparison utilities work", {
  data("toy_expr", package = "GLEAM")
  s1 <- score_pathway(expr = toy_expr$expr, meta = toy_expr$meta, geneset = "immune_small", seurat = FALSE, method = "rank", min_genes = 3, verbose = FALSE)
  s2 <- score_pathway(expr = toy_expr$expr, meta = toy_expr$meta, geneset = "immune_small", seurat = FALSE, method = "mean", min_genes = 3, verbose = FALSE)

  long <- pivot_scores_long(s1)
  expect_true(all(c("cell_id", "pathway", "score") %in% colnames(long)))

  cmp <- compare_scoring_methods(rank = s1, mean = s2)
  expect_true(is.data.frame(cmp))

  f <- tempfile(fileext = ".csv")
  export_scores(s1, f, format = "csv", include_meta = TRUE)
  expect_true(file.exists(f))
})
