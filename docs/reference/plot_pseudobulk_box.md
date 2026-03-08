# Pseudobulk boxplot for signature scores

Pseudobulk boxplot for signature scores

## Usage

``` r
plot_pseudobulk_box(
  score,
  pathway,
  group,
  sample,
  celltype = NULL,
  palette = "gleam_discrete",
  point_size = 1.6,
  alpha = 0.75,
  theme_params = list()
)
```

## Arguments

- score:

  `gleam_score` object.

- pathway:

  Signature name (legacy argument name).

- group:

  Group variable.

- sample:

  Sample variable.

- celltype:

  Optional celltype variable.

- palette:

  Discrete palette name or custom colors.

- point_size:

  Jitter point size.

- alpha:

  Box alpha.

- theme_params:

  Optional list passed to
  [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md).

## Value

A `ggplot` object.
