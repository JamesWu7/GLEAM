test_that("geneset parser handles built-in and list input", {
  gs <- get_geneset("immune_small")
  expect_true(is.list(gs))
  expect_true(length(gs) >= 1)

  custom <- list(A = c("G1", "G2"), B = c("G3", "G4"))
  gs2 <- as_geneset(custom)
  expect_equal(names(gs2), c("A", "B"))
})

test_that("geneset matching returns overlap info", {
  gs <- list(P1 = c("A", "B", "C"), P2 = c("D", "E"))
  m <- suppressWarnings(match_geneset(gs, expr_genes = c("A", "C", "D"), verbose = FALSE))
  expect_true(all(c("geneset", "info") %in% names(m)))
  expect_equal(nrow(m$info), 2)
})
