# Box plot of sample-level signature scores

Box plot of sample-level signature scores

## Usage

``` r
plot_box(
  score,
  signature,
  group,
  sample,
  palette = "gleam_discrete",
  point_size = 1.8,
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

- sample:

  Sample variable (metadata column name or vector).

- palette:

  Discrete palette name or custom colors.

- point_size:

  Jitter point size.

- alpha:

  Box alpha.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
