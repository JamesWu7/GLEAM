# AUC-like pathway scoring on matrix input

AUC-like pathway scoring on matrix input

## Usage

``` r
score_auc_matrix(expr, genesets, auc_max_rank = 0.05, verbose = TRUE)
```

## Arguments

- expr:

  Gene-by-cell matrix.

- genesets:

  Named list of pathways.

- auc_max_rank:

  Fraction of top-ranked genes considered active.

- verbose:

  Whether to print progress.

## Value

Pathway-by-cell numeric matrix.
