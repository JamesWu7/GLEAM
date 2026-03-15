# Canonical signature testing interface

Canonical signature testing interface

## Usage

``` r
test_signature(
  score,
  group = NULL,
  sample = NULL,
  celltype = NULL,
  target_celltype = NULL,
  region = NULL,
  pseudotime = NULL,
  lineage = NULL,
  level = c("cell", "sample", "celltype", "sample_celltype", "pseudobulk", "region",
    "sample_region", "trajectory"),
  method = c("wilcox", "t", "lm", "limma", "edgeR", "DESeq2", "spearman", "tradeSeq"),
  ref_group = NULL,
  paired = FALSE,
  covariates = NULL,
  adjust_method = "BH",
  min_cells = 10,
  min_samples = 2,
  aggregation = c("mean", "median", "fraction", "sum"),
  threshold = 0,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot"),
  verbose = TRUE
)
```

## Arguments

- score:

  `gleam_score` object.

- group:

  Group variable (column name or vector).

- sample:

  Sample variable (column name or vector).

- celltype:

  Celltype variable (column name or vector).

- target_celltype:

  Target celltype label for within-celltype group comparisons.

- region:

  Spatial region/domain variable (column name or vector).

- pseudotime:

  Pseudotime source for trajectory mode.

- lineage:

  Lineage source for trajectory mode.

- level:

  Comparison level.

- method:

  Statistical method.

- ref_group:

  Reference group.

- paired:

  Placeholder, currently ignored.

- covariates:

  Placeholder, currently ignored.

- adjust_method:

  Multiple testing adjustment method.

- min_cells:

  Minimum cells per group.

- min_samples:

  Minimum samples per group.

- aggregation:

  Aggregation summary used by pseudobulk levels.

- threshold:

  Threshold used when `aggregation = 'fraction'`.

- backend:

  Trajectory backend used when `level = 'trajectory'`.

- verbose:

  Print messages.

## Value

An object of class `gleam_test`.
