# GLEAM

## GLEAM: Gene-set and cell-state exploration across space and time in R

[![R-CMD-check](https://github.com/jameswoo/GLEAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jameswoo/GLEAM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/jameswoo/GLEAM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/jameswoo/GLEAM/actions/workflows/pkgdown.yaml)

- **G** = Gene-set
- **L** = cell-state
- **E** = Exploration
- **A** = across
- **M** = space and time

GLEAM is a matrix-first toolkit for pathway and cell-state analysis across single-cell and spatial transcriptomics.

> If your repository still uses the old working name, rename it to **`GLEAM`**.

## Installation

### Fresh R setup
```r
source("scripts/setup_dev.R")
source("scripts/check_env.R")
source("scripts/check_optional_deps.R")
```

### Build and test
```r
devtools::document()
devtools::test()
devtools::check()
devtools::build()
```

## Dependency strategy
- Core workflow imports only lightweight packages (`Matrix`, `ggplot2`, base R stats/utils).
- Rich integrations are optional (Seurat, msigdbr, trajectory packages, wrapper scorers, palette/plot extras).
- Missing optional dependencies fail with clear install guidance.

## Quickstart (matrix workflow)
```r
library(GLEAM)
data("pbmc_medium_matrix", package = "GLEAM")
data("pbmc_medium_meta", package = "GLEAM")

sc <- score_pathway(
  expr = pbmc_medium_matrix,
  meta = pbmc_medium_meta,
  geneset = "immune_small",
  geneset_source = "builtin",
  seurat = FALSE,
  method = "ensemble"
)
```

## Seurat v4 / v5 workflow
```r
# Requires Seurat / SeuratObject
# sc <- score_pathway(
#   object = seu,
#   geneset = "hallmark",
#   seurat = TRUE,
#   assay = "RNA",
#   layer = "data", # v5 layer path
#   slot = "data",  # v4 slot fallback
#   method = "rank"
# )
```

## Geneset workflow
```r
list_geneset_sources()
search_geneset("hallmark", "INTERFERON", source = "builtin")

# built-in / list / GMT / data.frame
gs_builtin <- get_geneset("hallmark", source = "builtin")

# optional msigdbr-backed sources
# gs_go <- get_geneset(NULL, source = "go", ontology = "BP")
# gs_kegg <- get_geneset(NULL, source = "kegg")
# gs_reactome <- get_geneset(NULL, source = "reactome")
```

## Multiple scoring methods
```r
list_scoring_methods()

sc_rank <- score_pathway(expr = pbmc_medium_matrix, meta = pbmc_medium_meta,
                         geneset = "immune_small", seurat = FALSE, method = "rank")
sc_auc <- score_pathway(expr = pbmc_medium_matrix, meta = pbmc_medium_meta,
                        geneset = "immune_small", seurat = FALSE, method = "auc")

compare_scoring_methods(rank = sc_rank, auc = sc_auc)
```

## Differential analysis

### Pseudobulk
```r
res_pb <- test_pathway(
  score = sc,
  group = "group",
  sample = "sample",
  celltype = "celltype",
  level = "pseudobulk",
  method = "wilcox"
)
```

### Within-celltype
```r
res_wct <- compare_groups_within_celltype(
  score = sc,
  group = "group",
  celltype = "celltype",
  target_celltype = "CD8_T",
  level = "sample"
)
```

## Trajectory mapping and differential
```r
res_traj <- test_pathway_trajectory(sc, pseudotime = "pseudotime", lineage = "lineage")
plot_pseudotime_score(sc, pathway = rownames(sc$score)[1], pseudotime = "pseudotime", lineage = "lineage")
```

## Spatial scoring and differential
```r
data("spatial_medium_expr", package = "GLEAM")
data("spatial_medium_meta", package = "GLEAM")
data("spatial_medium_coords", package = "GLEAM")

sc_sp <- score_pathway(expr = spatial_medium_expr, meta = spatial_medium_meta,
                       geneset = "immune_small", seurat = FALSE, method = "rank")
plot_spatial_score(sc_sp, pathway = rownames(sc_sp$score)[1], coords = spatial_medium_coords)

res_sp <- test_pathway_spatial(sc_sp, region = "region", group = "group", sample = "sample", level = "sample_region")
```

## Visualization examples
```r
plot_violin(sc, pathway = rownames(sc$score)[1], group = "group")
plot_box(sc, pathway = rownames(sc$score)[1], group = "group", sample = "sample")
plot_dot(sc, by = c("group", "celltype"))
plot_heatmap(sc, by = c("group", "celltype"), top_n = 20)
plot_volcano(res_pb)
```

## Export and reporting
```r
long_scores <- pivot_scores_long(sc)
export_scores(sc, file = "gleam_scores.csv", format = "csv", include_meta = TRUE)
```

## Tutorials and local HTML rendering
```r
source("scripts/render_examples.R")
```
Rendered files: `docs/tutorial_html/`.

## Documentation publishing
- pkgdown config: `_pkgdown.yml`
- published site target: https://jameswoo.github.io/GLEAM
