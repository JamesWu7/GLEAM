# Rank-based pathway scoring on matrix input

Rank-based pathway scoring on matrix input

## Usage

``` r
score_rank_matrix(expr, genesets, normalize = TRUE, verbose = TRUE)
```

## Arguments

- expr:

  Gene-by-cell matrix.

- genesets:

  Named list of pathways.

- normalize:

  Whether to normalize rank score to the 0-1 range.

- verbose:

  Whether to print progress.

## Value

Pathway-by-cell numeric matrix.
