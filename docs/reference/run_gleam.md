# End-to-end GLEAM workflow

End-to-end GLEAM workflow

## Usage

``` r
run_gleam(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset = "hallmark",
  geneset_source = "auto",
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = "ensemble",
  group = NULL,
  sample = NULL,
  celltype = NULL,
  comparison = c("group", "celltype", "within_celltype", "trajectory", "spatial"),
  target_celltype = NULL,
  level = "sample",
  pseudotime = NULL,
  lineage = NULL,
  region = NULL,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot"),
  top_n = 20,
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat object.

- expr:

  Expression matrix for matrix mode.

- meta:

  Metadata for matrix mode.

- geneset:

  Geneset input.

- geneset_source:

  Geneset source.

- seurat:

  Input mode.

- assay:

  Seurat assay.

- layer:

  Seurat layer.

- slot:

  Legacy Seurat slot fallback.

- method:

  Scoring method.

- group:

  Group variable.

- sample:

  Sample variable.

- celltype:

  Celltype variable.

- comparison:

  Comparison mode.

- target_celltype:

  Target celltype for within-celltype comparison.

- level:

  Testing level.

- pseudotime:

  Pseudotime vector or metadata column name.

- lineage:

  Lineage vector or metadata column name.

- region:

  Spatial region/domain column name or vector.

- backend:

  Trajectory backend used when `comparison = 'trajectory'`.

- top_n:

  Number of top pathways.

- verbose:

  Print messages.

## Value

List containing score object and comparison result.
