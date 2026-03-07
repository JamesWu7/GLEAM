test_that("matrix mode input checks work", {
  data("toy_expr", package = "GLEAM")
  expect_true(GLEAM:::check_input_mode(object = NULL, expr = toy_expr$expr, seurat = FALSE))
  expect_error(GLEAM:::check_input_mode(object = NULL, expr = NULL, seurat = FALSE), "expr")
})

test_that("seurat mode validates Seurat class", {
  fake <- list()
  class(fake) <- "NotSeurat"
  expect_error(GLEAM:::check_input_mode(object = fake, expr = NULL, seurat = TRUE), "Seurat")
})
