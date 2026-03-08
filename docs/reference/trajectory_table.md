# Build trajectory table from score object

Build trajectory table from score object

## Usage

``` r
trajectory_table(
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

Data.frame.
