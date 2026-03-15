# Split violin plot for signature scores

Split violin plot for signature scores

## Usage

``` r
plot_split_violin(
  score,
  signature = NULL,
  x,
  split.by,
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

- x:

  Grouping variable for x-axis.

- split.by:

  Variable used for split/fill.

- palette:

  Discrete palette name or custom colors.

- alpha:

  Violin transparency.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
