test_that("public API keeps canonical signature entrypoints", {
  ns <- getNamespaceExports("GLEAM")

  expect_true(all(c(
    "run_gleam",
    "score_signature",
    "aggregate_signature",
    "test_signature",
    "test_signature_spatial",
    "test_signature_trajectory",
    "compare_scoring_methods"
  ) %in% ns))
})

test_that("legacy and duplicate entrypoints are no longer exported", {
  ns <- getNamespaceExports("GLEAM")

  expect_false(any(c(
    "run_scpathway",
    "score_pathway",
    "aggregate_pathway",
    "test_pathway",
    "test_pathway_spatial",
    "test_pathway_trajectory",
    "differential_pathway",
    "differential_pathway_spatial",
    "differential_pathway_trajectory",
    "differential_signature",
    "differential_signature_spatial",
    "differential_signature_trajectory",
    "compare_methods"
  ) %in% ns))
})
