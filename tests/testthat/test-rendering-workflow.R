test_that("render script exists", {
  candidates <- c(
    file.path("scripts", "render_examples.R"),
    file.path("..", "scripts", "render_examples.R"),
    file.path("..", "..", "scripts", "render_examples.R")
  )
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f)) {
    skip("render_examples.R is not bundled in this check context")
  }
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_true(grepl("rmarkdown::render", txt, fixed = TRUE))
})
