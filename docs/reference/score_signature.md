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
  method = c("rank", "mean", "zscore", "scaled_mean", "robust_mean", "ensemble",
    "AddModuleScore", "UCell", "AUCell", "ssGSEA", "GSVA", "singscore"),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "zscore", "mean"),
  ensemble_combine = c("mean", "median"),
  ensemble_standardize = c("zscore", "rank"),
  ensemble_weights = NULL,
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

  Geneset input. Also supports signed signatures:
  `list(signatureA = list(up = c(...), down = c(...)))`.

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

  Scoring method. Canonical methods include: `rank`, `mean`, `zscore`,
  `scaled_mean`, `robust_mean`, `ensemble`, `AddModuleScore`, `UCell`,
  `AUCell`, `ssGSEA`, `GSVA`, `singscore`.

- min_genes:

  Minimum genes per pathway.

- max_genes:

  Maximum genes per pathway.

- auc_max_rank:

  AUC top-rank proportion used by `AUCell`.

- ensemble_methods:

  Methods used by ensemble.

- ensemble_combine:

  Ensemble combine strategy.

- ensemble_standardize:

  Harmonization before ensemble aggregation: `zscore` (recommended) or
  `rank`.

- ensemble_weights:

  Optional named numeric vector of per-method weights.

- method_params:

  Optional named list for method-specific parameters. Supported keys
  include: `auc_max_rank` (AUCell), `alpha` or `ssgsea_alpha` (ssGSEA),
  `nbin`/`ctrl`/`seed` (AddModuleScore), `kcdf` (GSVA), `ucell_max_rank`
  (UCell).

- verbose:

  Whether to print messages.

## Value

An object of class `gleam_score`.

## Method Guide

|                  |                                |                                                             |
|------------------|--------------------------------|-------------------------------------------------------------|
| Method           | Best for                       | Notes                                                       |
| `rank`           | Fast robust baseline           | Rank aggregation concept aligned with singscore             |
| `UCell`          | Large single-cell data         | Mann-Whitney U-based robust scoring                         |
| `AUCell`         | Rank-enrichment scoring        | SCENIC/AUCell implementation                                |
| `AddModuleScore` | Seurat-native workflows        | Signature minus control-bin expression                      |
| `ssGSEA`         | Bulk-like enrichment workflows | GSVA package ssGSEA implementation                          |
| `GSVA`           | Bulk expression matrices       | GSVA kernel-based pathway variation                         |
| `mean`           | Simple signatures              | Mean expression aggregation                                 |
| `zscore`         | Cross-signature comparability  | Per-gene z-score then aggregate                             |
| `robust_mean`    | Outlier-prone data             | Median/MAD-scaled robust aggregation                        |
| `ensemble`       | Consensus scoring              | Harmonize each method (`zscore` or `rank`) before averaging |

## Method Parameters (`method_params`)

Supported keys include: `auc_max_rank` (AUCell), `alpha` or
`ssgsea_alpha` (ssGSEA), `nbin`/`ctrl`/`seed` (AddModuleScore), `kcdf`
(GSVA), `ucell_max_rank` (UCell). For ensemble, use `ensemble_methods`,
`ensemble_standardize`, and optional `ensemble_weights`.
