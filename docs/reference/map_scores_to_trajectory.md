# Map scores to trajectory table

Map scores to trajectory table

## Usage

``` r
map_scores_to_trajectory(
  score,
  pseudotime = NULL,
  lineage = NULL,
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

- backend:

  Trajectory backend. `auto` detects from provided inputs.

## Value

Long-format data.frame.
