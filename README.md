<p align="center">
  <img src="man/figures/GLEAM_LOG.jpg" alt="GLEAM logo" width="220"/>
</p>

# GLEAM

GLEAM: Gene-set and cell-state exploration across space and time in R

[![R-CMD-check](https://github.com/jameswoo/GLEAM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jameswoo/GLEAM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/jameswoo/GLEAM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/jameswoo/GLEAM/actions/workflows/pkgdown.yaml)

**Navigation:** [Documentation](https://jameswoo.github.io/GLEAM/) | [Reference](https://jameswoo.github.io/GLEAM/reference/) | [Tutorials](https://jameswoo.github.io/GLEAM/articles/) | [Citation](https://jameswoo.github.io/GLEAM/articles/GLEAM_citation.html)

GLEAM provides pathway and cell-state analysis workflows across scRNA-seq, spatial transcriptomics, and pseudotime/trajectory settings. It supports matrix-native and Seurat-based analysis, unified pathway scoring, differential testing after scoring, trajectory/spatial mapping, publication-ready plots, and tutorial-first onboarding. Built-in workflows focus on **human** and **mouse** gene sets, while custom gene sets remain open for other species.

## Quick start (Seurat scRNA-seq)

```r
library(GLEAM)

# seu: preprocessed Seurat object (v4/v5)
hallmark_gs <- get_geneset("hallmark", source = "builtin", species = "human")
kegg_gs <- get_geneset(NULL, source = "kegg", species = "human")  # optional msigdbr backend

sc_h <- score_pathway(
  object = seu,
  geneset = hallmark_gs,
  geneset_source = "list",
  seurat = TRUE,
  assay = "RNA",
  layer = "data",
  slot = "data",
  method = "ensemble"
)

res_h <- test_pathway(
  score = sc_h,
  group = "group",
  sample = "sample",
  celltype = "celltype",
  level = "pseudobulk",
  method = "wilcox"
)

# Seurat-style practical displays
top_pw <- res_h$table$pathway[order(res_h$table$p_adj)][1]
plot_embedding_score(sc_h, pathway = top_pw, object = seu, reduction = "umap")
plot_embedding_score(sc_h, pathway = top_pw, object = seu, reduction = "pca")
plot_embedding_score(sc_h, pathway = top_pw, object = seu, reduction = "tsne")
plot_dot_bar(sc_h, by = c("group", "celltype"))
plot_violin(sc_h, pathway = rownames(sc_h$score)[1], group = "group")
plot_ridge(sc_h, pathway = rownames(sc_h$score)[1], group = "celltype")
plot_volcano(res_h)

# Optional KEGG run
sc_k <- score_pathway(object = seu, geneset = kegg_gs, geneset_source = "list", seurat = TRUE, method = "rank")
```

## Quick start (Seurat spatial transcriptomics)

```r
library(GLEAM)

# sp_seu: Seurat spatial object with coordinates and image when available
gs <- get_geneset("hallmark", source = "builtin", species = "human")

sp <- score_pathway(
  object = sp_seu,
  geneset = gs,
  geneset_source = "list",
  seurat = TRUE,
  assay = "Spatial",
  layer = "data",
  slot = "data",
  method = "rank"
)

# Slice/tissue view when image exists
coords <- data.frame(
  x = sp_seu@meta.data$x,
  y = sp_seu@meta.data$y,
  row.names = rownames(sp_seu@meta.data)
)
img <- as.raster(matrix(colorRampPalette(c("#f7f3e8", "#eadfca", "#d9c7a4"))(256), nrow = 16))
plot_spatial_score(sp, pathway = rownames(sp$score)[1], coords = coords, image = img, split.by = "sample")

sp_res <- test_pathway(
  score = sp,
  region = "region",
  group = "group",
  sample = "sample",
  level = "sample_region",
  method = "wilcox"
)
top_sp_pw <- sp_res$table$pathway[order(sp_res$table$p_adj)][1]
plot_spatial_score(sp, pathway = top_sp_pw, coords = coords, image = img, split.by = "region")
plot_spatial_compare(sp_res)
```

## Custom gene-set example (concise)

```r
custom_gs <- list(
  IFN_custom = c("STAT1", "IRF1", "ISG15", "IFIT3"),
  CYT_custom = c("NKG7", "PRF1", "GZMB", "GNLY")
)

sc_custom <- score_pathway(object = seu, geneset = custom_gs, geneset_source = "list", seurat = TRUE, method = "mean")
plot_dot(sc_custom, by = c("group", "celltype"))
```

## Supported gene-set sources

- `builtin` (human/mouse focus): Hallmark-like and immune sets shipped in package.
- `list`: named list input.
- `gmt`: GMT file parser.
- `data.frame`: `pathway` + `gene` columns.
- `msigdb`, `go`, `kegg`, `reactome`: optional `msigdbr` backend, explicitly dependency-gated.

## Visualization parameters

All plotting functions return `ggplot` and support shared readability controls through `gleam_theme()`:

- Text/font: `base_size`, `title_size`, `axis_text_size`, `legend_text_size`, `font_family`, `font_face`.
- Text colors: `title_color`, `axis_text_color`, `legend_title_color`, `text_color`.
- Palette/color: `palette` arguments in plot functions plus `get_palette()`, `scale_gleam_color()`, `scale_gleam_fill()`.
- Common chart controls: `point_size`, `alpha`, `split.by`, `reduction`, and spatial image vs coordinate-only displays.

Example:

```r
p <- plot_violin(sc_h, pathway = rownames(sc_h$score)[1], group = "group")
apply_gleam_theme(p, base_size = 14, title_size = 18, axis_text_size = 12, font_family = "sans")
```

## Detailed tutorials by function category

- Input and extraction: `score_pathway()`, `extract_embedding()`, `seurat_mode()`, `spatial_table()`.
- Genesets: `get_geneset()`, `list_geneset_sources()`, `search_geneset()`, `as_geneset()`, `read_gmt()`.
- Scoring: `score_pathway()`, `list_scoring_methods()`, `compare_scoring_methods()`.
- Differential analysis: `test_pathway()`, `compare_celltypes()`, `compare_groups_within_celltype()`.
- Trajectory: `test_pathway_trajectory()`, `plot_pseudotime_score()`, `plot_trajectory_score()`.
- Spatial: `test_pathway_spatial()`, `plot_spatial_score()`, `plot_spatial_compare()`, `plot_spatial_multi()`.
- Visualization: `plot_dot_bar()`, `plot_violin()`, `plot_ridge()`, `plot_embedding_score()`, `plot_volcano()`.
- Export/comparison: `collect_scores()`, `summarize_scores()`, `pivot_scores_long()`, `export_scores()`.

## Full workflow tutorials

- Full scRNA-seq workflow: [GLEAM_full_scrna_workflow](https://jameswoo.github.io/GLEAM/articles/GLEAM_full_scrna_workflow.html)
- Full spatial workflow: [GLEAM_full_spatial_workflow](https://jameswoo.github.io/GLEAM/articles/GLEAM_full_spatial_workflow.html)

## Citation

- R-native citation: `citation("GLEAM")`
- Citation page: [GLEAM citation guide](https://jameswoo.github.io/GLEAM/articles/GLEAM_citation.html)

## Local tutorial rendering

```r
source("scripts/render_examples.R")
```

Rendered HTML output directory: `docs/tutorial_html/`.
