# Plot multiple signatures on spatial coordinates

Plot multiple signatures on spatial coordinates

## Usage

``` r
plot_spatial_multi(
  score,
  pathways,
  coords = NULL,
  object = NULL,
  palette = "gleam_continuous",
  point_size = 1.2,
  alpha = 0.9,
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- pathways:

  Character vector of signature names.

- coords:

  Spatial coordinates data.frame with x/y.

- object:

  Optional Seurat spatial object for in-slice plotting.

- palette:

  Continuous palette.

- point_size:

  Point size.

- alpha:

  Point alpha.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
