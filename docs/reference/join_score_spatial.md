# Join score with spatial coordinates

Join score with spatial coordinates

## Usage

``` r
join_score_spatial(score, coords, meta = NULL)
```

## Arguments

- score:

  `gleam_score` object.

- coords:

  Spatial coordinates data.frame or matrix with x/y.

- meta:

  Optional metadata to merge.

## Value

Data.frame containing cell_id, coordinates, and metadata.
