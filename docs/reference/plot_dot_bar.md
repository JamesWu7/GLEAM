# Dot + bar summary plot for signature groups

Dot + bar summary plot for signature groups

## Usage

``` r
plot_dot_bar(
  score,
  by,
  threshold = 0,
  pathway = NULL,
  color_palette = "gleam_continuous",
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

- pathway:

  Optional signature filter (legacy argument name).

- color_palette:

  Continuous palette for mean score.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
