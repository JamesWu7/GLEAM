# Coerce geneset input to named list

Coerce geneset input to named list

## Usage

``` r
as_geneset(x)
```

## Arguments

- x:

  Geneset input.

## Value

Named list. Supports standard pathways (`list(pathway = c(...))`) and
signed pathways (`list(pathway = list(up = c(...), down = c(...)))`).
