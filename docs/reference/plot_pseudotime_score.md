# Plot signature score over pseudotime

Plot signature score over pseudotime

## Usage

``` r
plot_pseudotime_score(
  score,
  signature = NULL,
  pseudotime = NULL,
  lineage = NULL,
  smooth = TRUE,
  point_size = 1.1,
  alpha = 0.55,
  palette = "gleam_discrete",
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- signature:

  Signature name.

- pseudotime:

  Pseudotime source.

- lineage:

  Optional lineage source for coloring.

- smooth:

  Add smoothing line.

- point_size:

  Point size.

- alpha:

  Point alpha.

- palette:

  Discrete palette for lineages.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://jameswu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
