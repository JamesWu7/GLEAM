# Join score matrix and trajectory metadata

Join score matrix and trajectory metadata

## Usage

``` r
join_score_trajectory(
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
