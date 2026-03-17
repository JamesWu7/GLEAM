# Extract spatial coordinates

Extract spatial coordinates

## Usage

``` r
extract_spatial_coords(
  object = NULL,
  meta = NULL,
  coords = NULL,
  seurat = TRUE
)
```

## Arguments

- object:

  Seurat object.

- meta:

  Metadata for matrix mode.

- coords:

  Optional direct coordinate data.frame.

- seurat:

  Input mode.

## Value

Data.frame with columns `x` and `y`.
