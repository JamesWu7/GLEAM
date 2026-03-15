# Canonical trajectory signature testing interface

Canonical trajectory signature testing interface

## Usage

``` r
test_signature_trajectory(
  score,
  pathway = NULL,
  pseudotime = NULL,
  lineage = NULL,
  method = c("spearman", "lm", "tradeSeq"),
  adjust_method = "BH",
  verbose = TRUE,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
)
```

## Arguments

- score:

  `gleam_score` object.

- pathway:

  Optional pathway name. `NULL` tests all pathways.

- pseudotime:

  Pseudotime source.

- lineage:

  Lineage source.

- method:

  One of `spearman`, `lm`, `tradeSeq`.

- adjust_method:

  P-value adjustment method.

- verbose:

  Print messages.

- backend:

  Trajectory backend to use. `auto` detects from provided inputs.

## Value

`gleam_test` object.
