# Plot signature score on spatial coordinates

Plot signature score on spatial coordinates

## Usage

``` r
plot_spatial_score(
  score,
  signature,
  coords = NULL,
  object = NULL,
  image = NULL,
  palette = "gleam_continuous",
  split.by = NULL,
  point_size = 1.4,
  alpha = 0.9,
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- signature:

  Signature name.

- coords:

  Spatial coordinates data.frame with x/y.

- object:

  Optional Seurat spatial object for in-slice plotting.

- image:

  Optional raster/image object. If provided, used as background.

- palette:

  Continuous palette name or colors.

- split.by:

  Optional metadata split variable.

- point_size:

  Point size.

- alpha:

  Point alpha.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
