# Canonical signature scoring interface

Canonical signature scoring interface

## Usage

``` r
score_signature(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset,
  geneset_source = c("auto", "builtin", "list", "gmt", "data.frame", "msigdb", "go",
    "kegg", "reactome"),
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  ontology = c("BP", "MF", "CC"),
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = c("rank", "auc", "zmean", "mean", "scaled_mean", "robust", "singscore_like",
    "ssgsea_like", "ensemble", "addmodulescore", "ucell", "aucell", "gsva", "singscore"),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "auc", "zmean"),
  ensemble_combine = c("mean", "median"),
  method_params = list(),
  verbose = TRUE
)
```

## Arguments

- object:

  Seurat object when `seurat = TRUE`.

- expr:

  Expression matrix or `dgCMatrix` when `seurat = FALSE`.

- meta:

  Metadata data.frame for matrix mode.

- geneset:

  Geneset input.

- geneset_source:

  Geneset source (`auto`, `builtin`, `list`, `gmt`, `data.frame`,
  `msigdb`, `go`, `kegg`, `reactome`).

- species:

  Species label for source-aware geneset loading.

- collection:

  Collection parameter for source-aware geneset loading.

- subcollection:

  Subcollection parameter for source-aware geneset loading.

- ontology:

  GO ontology (`BP`, `MF`, `CC`) for GO source.

- seurat:

  Logical input mode.

- assay:

  Seurat assay name.

- layer:

  Seurat v5 layer name. If `NULL`, tries `data` then `counts`.

- slot:

  Legacy Seurat slot fallback.

- method:

  Scoring method. See
  [`list_scoring_methods()`](https://JamesWu7.github.io/GLEAM/reference/list_scoring_methods.md).

- min_genes:

  Minimum genes per pathway.

- max_genes:

  Maximum genes per pathway.

- auc_max_rank:

  AUC top-rank proportion.

- ensemble_methods:

  Methods used by ensemble.

- ensemble_combine:

  Ensemble combine strategy.

- method_params:

  Optional named list for method-specific parameters.

- verbose:

  Whether to print messages.

## Value

An object of class `gleam_score`.
