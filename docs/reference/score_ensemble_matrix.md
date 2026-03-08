# Ensemble pathway scoring on matrix input

Ensemble pathway scoring on matrix input

## Usage

``` r
score_ensemble_matrix(
  expr,
  genesets,
  methods = c("rank", "auc", "zmean"),
  combine = c("mean", "median"),
  auc_max_rank = 0.05,
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

- auc_max_rank:

  Fraction of top ranked genes for auc method.

- verbose:

  Whether to print progress.

## Value

Pathway-by-cell numeric matrix.
