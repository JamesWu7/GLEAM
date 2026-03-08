# Heatmap of aggregated signature scores

Heatmap of aggregated signature scores

## Usage

``` r
plot_heatmap(
  score,
  by,
  fun = c("mean", "median", "fraction"),
  threshold = 0,
  top_n = NULL,
  palette = "gleam_continuous",
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- by:

  Grouping metadata columns.

- fun:

  Aggregation function.

- threshold:

  Threshold used by `fraction`.

- top_n:

  Optional top signatures by variance.

- palette:

  Continuous palette name or custom colors.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
