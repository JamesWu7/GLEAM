test_that("additional plotting APIs return ggplot", {
  data("toy_expr", package = "GLEAM")
  sc <- score_signature(expr = toy_expr$expr, meta = toy_expr$meta, geneset = "immune_small", seurat = FALSE, method = "rank", min_genes = 3, verbose = FALSE)

  p1 <- plot_split_violin(sc, pathway = rownames(sc$score)[1], x = "celltype", split.by = "group")
  expect_s3_class(p1, "ggplot")

  p2 <- plot_dot_bar(sc, by = c("group", "celltype"), pathway = rownames(sc$score)[1])
  expect_s3_class(p2, "ggplot")

  sc2 <- sc
  sc2$meta$stim <- sc2$meta$group
  sc2$meta$seurat_annotations <- sc2$meta$celltype
  p2b <- plot_dot_bar(sc2, by = c("stim", "seurat_annotations"), pathway = rownames(sc2$score)[1])
  expect_s3_class(p2b, "ggplot")

  p3 <- plot_pseudobulk_box(sc, pathway = rownames(sc$score)[1], group = "group", sample = "sample", celltype = "celltype")
  expect_s3_class(p3, "ggplot")
})

test_that("palette helpers are available", {
  pals <- list_palettes()
  expect_true("gleam_discrete" %in% pals)
  cols <- get_palette("gleam_discrete", n = 5)
  expect_equal(length(cols), 5)
})

test_that("theme helper returns ggplot theme", {
  th <- gleam_theme(base_size = 12, title_size = 16, axis_text_size = 11)
  expect_s3_class(th, "theme")
})
