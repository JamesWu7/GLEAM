# Volcano plot for differential signatures

Volcano plot for differential signatures

## Usage

``` r
plot_volcano(
  x,
  p_col = "p_adj",
  effect_col = "effect_size",
  sig_thresh = 0.05,
  point_size = 1.8,
  alpha = 0.8,
  theme_params = list()
)
```

## Arguments

- x:

  `gleam_test` object or result data.frame.

- p_col:

  P-value column name.

- effect_col:

  Effect size column name.

- sig_thresh:

  Significance threshold for adjusted p-value.

- point_size:

  Point size.

- alpha:

  Point alpha.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
