# ssGSEA-like pathway scoring on matrix input

ssGSEA-like pathway scoring on matrix input

## Usage

``` r
score_ssgsea_like_matrix(expr, genesets, alpha = 0.25, verbose = TRUE)
```

## Arguments

- expr:

  Gene-by-cell matrix.

- genesets:

  Named list of pathways.

- alpha:

  Rank-weight exponent.

- verbose:

  Whether to print progress.

## Value

Pathway-by-cell numeric matrix.
