# scPathway

`scPathway` is a lightweight, low-dependency package for pathway analysis in single-cell RNA-seq data.

## Features
- Single-cell pathway scoring (`rank`, `auc`, `zmean`, `ensemble`)
- Cell-level exploratory comparison
- Sample-level and sample+celltype formal comparison
- Celltype comparison and within-celltype group comparison
- Seurat v5 optional compatibility and matrix-native workflow

## First-run setup (fresh R)
```r
source("scripts/setup_dev.R")
source("scripts/check_env.R")
```

## Dev commands
```r
devtools::document()
devtools::test()
devtools::check()
devtools::build()
```

## Matrix workflow example
```r
library(scPathway)
data("toy_expr", package = "scPathway")

sc <- score_pathway(
  expr = toy_expr$expr,
  meta = toy_expr$meta,
  geneset = "immune_small",
  seurat = FALSE,
  method = "rank"
)

cmp <- compare_groups_within_celltype(
  score = sc,
  group = "group",
  celltype = "celltype",
  target_celltype = "CD8_T",
  level = "cell"
)

p <- plot_violin(sc, pathway = rownames(sc$score)[1], group = "group")
print(p)
```

## Seurat v5 workflow example
```r
# Requires Seurat/SeuratObject installed
# sc <- score_pathway(
#   object = seu,
#   geneset = "hallmark",
#   seurat = TRUE,
#   assay = "RNA",
#   layer = "data",
#   method = "ensemble"
# )
```
