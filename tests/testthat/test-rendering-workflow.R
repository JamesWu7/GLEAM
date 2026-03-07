test_that("render script exists", {
  f <- file.path("..", "..", "scripts", "render_examples.R")
  expect_true(file.exists(f))
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_true(grepl("rmarkdown::render", txt, fixed = TRUE))
})
