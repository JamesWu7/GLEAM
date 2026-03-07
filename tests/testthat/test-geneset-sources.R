test_that("geneset source listing works", {
  src <- list_geneset_sources()
  expect_true(is.data.frame(src))
  expect_true(all(c("source", "requires_package") %in% colnames(src)))
})

test_that("builtin geneset and search works", {
  gs <- get_geneset("hallmark", source = "builtin")
  expect_true(is.list(gs))
  hit <- search_geneset("hallmark", "INTERFERON", source = "builtin")
  expect_true(length(hit) >= 1)
})
