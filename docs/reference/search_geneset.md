# Search pathways in a geneset collection

Search pathways in a geneset collection

## Usage

``` r
search_geneset(geneset = "hallmark", pattern, ignore.case = TRUE, ...)
```

## Arguments

- geneset:

  Geneset input accepted by
  [`get_geneset()`](https://jameswu7.github.io/GLEAM/reference/get_geneset.md).

- pattern:

  Search pattern.

- ignore.case:

  Whether search is case-insensitive.

- ...:

  Additional arguments passed to
  [`get_geneset()`](https://jameswu7.github.io/GLEAM/reference/get_geneset.md).

## Value

Character vector of matched pathway names.
