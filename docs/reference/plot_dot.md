# Dot plot of signature summaries

Dot size represents fraction of cells above threshold, and dot color
represents mean signature score.

## Usage

``` r
plot_dot(
  score,
  by,
  threshold = 0,
  palette = "gleam_continuous",
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- by:

  Grouping metadata columns.

- threshold:

  Threshold for active fraction.

- palette:

  Continuous palette name or custom colors.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
