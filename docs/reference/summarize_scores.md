# Summarize score table

Summarize score table

## Usage

``` r
summarize_scores(
  score,
  by = c("sample", "group"),
  fun = c("mean", "median", "fraction"),
  threshold = 0
)
```

## Arguments

- score:

  `gleam_score` object.

- by:

  Grouping columns.

- fun:

  Aggregation function.

- threshold:

  Threshold used for fraction.

## Value

Summarized data.frame.
