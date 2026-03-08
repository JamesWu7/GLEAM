test_that("detect_seurat_version is robust", {
  expect_true(is.na(GLEAM:::detect_seurat_version(NULL)))

  fake <- structure(list(), class = "Seurat")
  v <- GLEAM:::detect_seurat_version(fake)
  expect_true(is.na(v) || is.numeric(v))
})

test_that("matrix mode works without Seurat", {
  data("toy_expr", package = "GLEAM")
  expect_true(GLEAM:::check_input_mode(object = NULL, expr = toy_expr$expr, seurat = FALSE))
})
