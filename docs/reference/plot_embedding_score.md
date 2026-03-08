# Plot signature score on embedding coordinates

Plot signature score on embedding coordinates

## Usage

``` r
plot_embedding_score(
  score,
  pathway,
  embedding = NULL,
  object = NULL,
  reduction = "umap",
  split.by = NULL,
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

- embedding:

  Embedding matrix with at least 2 columns.

- object:

  Optional Seurat object.

- reduction:

  Reduction name used when `embedding` is NULL.

- split.by:

  Optional metadata column for faceting.

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
