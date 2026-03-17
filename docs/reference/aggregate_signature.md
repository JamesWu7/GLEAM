# Canonical signature aggregation interface

Canonical signature aggregation interface

## Usage

``` r
aggregate_signature(
  score,
  by,
  fun = c("mean", "median", "fraction", "sum"),
  threshold = 0,
  long = TRUE
)
```

## Arguments

- score:

  `gleam_score` object.

- by:

  Grouping variable(s): metadata column name(s), or vector.

- fun:

  Aggregation function: mean, median, fraction, sum, or a function.

- threshold:

  Threshold used by `fraction`.

- long:

  Return long-format table if TRUE.

## Value

Aggregated table (long) or signature-by-group matrix (wide).
