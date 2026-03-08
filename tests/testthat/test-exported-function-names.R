test_that("canonical signature APIs are exported", {
  ns <- getNamespaceExports("GLEAM")
  expect_true(all(c(
    "score_signature",
    "aggregate_signature",
    "test_signature",
    "differential_signature",
    "test_signature_spatial",
    "differential_signature_spatial",
    "test_signature_trajectory",
    "differential_signature_trajectory"
  ) %in% ns))
})
