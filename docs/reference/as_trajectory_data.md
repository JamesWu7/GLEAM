# Convert trajectory inputs to a unified data.frame

Convert trajectory inputs to a unified data.frame

## Usage

``` r
as_trajectory_data(
  score,
  pseudotime = NULL,
  lineage = NULL,
  embeddings = NULL,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
)
```

## Arguments

- score:

  `gleam_score` object.

- pseudotime:

  Pseudotime source.

- lineage:

  Lineage source.

- embeddings:

  Optional embedding matrix.

- backend:

  Trajectory backend. `auto` detects from provided inputs.

## Value

Data.frame with `cell_id`, `pseudotime`, `lineage`, and embedding
columns.
