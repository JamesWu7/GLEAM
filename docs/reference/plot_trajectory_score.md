# Plot signature score on trajectory embedding

Plot signature score on trajectory embedding

## Usage

``` r
plot_trajectory_score(
  score,
  pathway,
  embeddings = NULL,
  reduction = "umap",
  object = NULL,
  point_size = 1.1,
  alpha = 0.9,
  palette = "gleam_continuous",
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- pathway:

  Signature name (legacy argument name).

- embeddings:

  Embedding matrix with at least 2 columns.

- reduction:

  Reduction name for Seurat extraction when `embeddings` is NULL.

- object:

  Optional Seurat object for reduction extraction.

- point_size:

  Point size.

- alpha:

  Point alpha.

- palette:

  Continuous palette name or custom colors.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
