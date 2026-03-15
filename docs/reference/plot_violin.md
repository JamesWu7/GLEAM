# Violin plot of signature scores

Violin plot of signature scores

## Usage

``` r
plot_violin(
  score,
  signature = NULL,
  group,
  celltype = NULL,
  trim = TRUE,
  palette = "gleam_discrete",
  point_size = 1.4,
  alpha = 0.7,
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- signature:

  Signature name.

- group:

  Group variable (metadata column name or vector).

- celltype:

  Optional celltype filter (single label or vector).

- trim:

  Passed to `geom_violin()`.

- palette:

  Discrete palette name or custom colors.

- point_size:

  Point size for boxplot outlier overlay.

- alpha:

  Violin transparency.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
