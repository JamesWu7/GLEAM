# Ensemble pathway scoring on matrix input

Ensemble pathway scoring on matrix input

## Usage

``` r
score_ensemble_matrix(
  expr,
  genesets,
  methods = c("rank", "zscore", "mean"),
  combine = c("mean", "median"),
  standardize = c("zscore", "rank"),
  weights = NULL,
  verbose = TRUE
)
```

## Arguments

- expr:

  Gene-by-cell matrix.

- genesets:

  Named list of pathways.

- methods:

  Methods to combine.

- combine:

  Combination strategy: mean or median.

- standardize:

  Ensemble harmonization (`zscore` or `rank`).

- weights:

  Optional named numeric vector of per-method weights.

- verbose:

  Whether to print progress.

## Value

Pathway-by-cell numeric matrix.
