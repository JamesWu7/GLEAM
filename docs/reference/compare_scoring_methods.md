# Compare multiple scoring methods

Compare multiple scoring methods

## Usage

``` r
compare_scoring_methods(..., pathway = NULL, summary_fun = c("mean", "median"))
```

## Arguments

- ...:

  Named score objects from different methods.

- pathway:

  Optional pathway subset.

- summary_fun:

  Summary function (`mean` or `median`).

## Value

Data.frame of per-method pathway summaries.
