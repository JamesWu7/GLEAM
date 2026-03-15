# Ridge plot for signature score distributions

Ridge plot for signature score distributions

## Usage

``` r
plot_ridge(
  score,
  signature = NULL,
  group,
  palette = "gleam_discrete",
  alpha = 0.75,
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- signature:

  Signature name.

- group:

  Group variable.

- palette:

  Discrete palette name or custom colors.

- alpha:

  Transparency of ridges.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
