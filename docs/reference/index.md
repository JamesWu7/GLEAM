# Package index

## Input and data extraction

- [`score_signature()`](https://JamesWu7.github.io/GLEAM/reference/score_signature.md)
  : Canonical signature scoring interface
- [`seurat_mode()`](https://JamesWu7.github.io/GLEAM/reference/seurat_mode.md)
  : Check Seurat compatibility mode
- [`extract_embedding()`](https://JamesWu7.github.io/GLEAM/reference/extract_embedding.md)
  : Extract embedding matrix

## Geneset management

- [`get_geneset()`](https://JamesWu7.github.io/GLEAM/reference/get_geneset.md)
  : Get geneset collection
- [`list_geneset_sources()`](https://JamesWu7.github.io/GLEAM/reference/list_geneset_sources.md)
  : List geneset sources
- [`search_geneset()`](https://JamesWu7.github.io/GLEAM/reference/search_geneset.md)
  : Search pathways in a geneset collection
- [`as_geneset()`](https://JamesWu7.github.io/GLEAM/reference/as_geneset.md)
  : Coerce geneset input to named list
- [`read_gmt()`](https://JamesWu7.github.io/GLEAM/reference/read_gmt.md)
  : Read GMT file
- [`check_geneset()`](https://JamesWu7.github.io/GLEAM/reference/check_geneset.md)
  : Validate genesets by size
- [`match_geneset()`](https://JamesWu7.github.io/GLEAM/reference/match_geneset.md)
  : Match genesets to expression genes

## Scoring

- [`aggregate_signature()`](https://JamesWu7.github.io/GLEAM/reference/aggregate_signature.md)
  : Canonical signature aggregation interface
- [`list_scoring_methods()`](https://JamesWu7.github.io/GLEAM/reference/list_scoring_methods.md)
  : List available scoring methods
- [`compare_scoring_methods()`](https://JamesWu7.github.io/GLEAM/reference/compare_scoring_methods.md)
  : Compare multiple scoring methods
- [`run_gleam()`](https://JamesWu7.github.io/GLEAM/reference/run_gleam.md)
  : End-to-end GLEAM workflow

## Differential analysis

- [`test_signature()`](https://JamesWu7.github.io/GLEAM/reference/test_signature.md)
  : Canonical signature testing interface
- [`test_signature_spatial()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_spatial.md)
  : Canonical spatial signature testing interface
- [`test_signature_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_trajectory.md)
  : Canonical trajectory signature testing interface

## Trajectory analysis

- [`extract_pseudotime()`](https://JamesWu7.github.io/GLEAM/reference/extract_pseudotime.md)
  : Extract pseudotime values
- [`extract_lineage()`](https://JamesWu7.github.io/GLEAM/reference/extract_lineage.md)
  : Extract lineage labels
- [`as_trajectory_data()`](https://JamesWu7.github.io/GLEAM/reference/as_trajectory_data.md)
  : Convert trajectory inputs to a unified data.frame
- [`join_score_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/join_score_trajectory.md)
  : Join score matrix and trajectory metadata
- [`map_scores_to_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/map_scores_to_trajectory.md)
  : Map scores to trajectory table
- [`plot_pseudotime_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_pseudotime_score.md)
  : Plot signature score over pseudotime
- [`plot_trajectory_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_trajectory_score.md)
  : Plot signature score on trajectory embedding

## Spatial analysis

- [`join_score_spatial()`](https://JamesWu7.github.io/GLEAM/reference/join_score_spatial.md)
  : Join score with spatial coordinates
- [`plot_spatial_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_score.md)
  : Plot signature score on spatial coordinates
- [`plot_spatial_compare()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_compare.md)
  : Plot spatial comparison result
- [`plot_spatial_multi()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_multi.md)
  : Plot multiple signatures on spatial coordinates

## Visualization

- [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md)
  : GLEAM plotting theme
- [`apply_gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/apply_gleam_theme.md)
  : Apply GLEAM theme to a ggplot
- [`plot_dot()`](https://JamesWu7.github.io/GLEAM/reference/plot_dot.md)
  : Dot plot of signature summaries
- [`plot_dot_bar()`](https://JamesWu7.github.io/GLEAM/reference/plot_dot_bar.md)
  : Dot + bar summary plot for signature groups
- [`plot_violin()`](https://JamesWu7.github.io/GLEAM/reference/plot_violin.md)
  : Violin plot of signature scores
- [`plot_split_violin()`](https://JamesWu7.github.io/GLEAM/reference/plot_split_violin.md)
  : Split violin plot for signature scores
- [`plot_ridge()`](https://JamesWu7.github.io/GLEAM/reference/plot_ridge.md)
  : Ridge plot for signature score distributions
- [`plot_box()`](https://JamesWu7.github.io/GLEAM/reference/plot_box.md)
  : Box plot of sample-level signature scores
- [`plot_pseudobulk_box()`](https://JamesWu7.github.io/GLEAM/reference/plot_pseudobulk_box.md)
  : Pseudobulk boxplot for signature scores
- [`plot_embedding_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_embedding_score.md)
  : Plot signature score on embedding coordinates
- [`plot_heatmap()`](https://JamesWu7.github.io/GLEAM/reference/plot_heatmap.md)
  : Heatmap of aggregated signature scores
- [`plot_volcano()`](https://JamesWu7.github.io/GLEAM/reference/plot_volcano.md)
  : Volcano plot for differential signatures
- [`plot_celltype_compare()`](https://JamesWu7.github.io/GLEAM/reference/plot_celltype_compare.md)
  : Plot celltype comparison result
- [`plot_group_in_celltype()`](https://JamesWu7.github.io/GLEAM/reference/plot_group_in_celltype.md)
  : Plot group comparison within a fixed celltype
- [`get_palette()`](https://JamesWu7.github.io/GLEAM/reference/get_palette.md)
  : Get palette colors
- [`list_palettes()`](https://JamesWu7.github.io/GLEAM/reference/list_palettes.md)
  : List available GLEAM palettes
- [`scale_gleam_color()`](https://JamesWu7.github.io/GLEAM/reference/scale_gleam_color.md)
  : GLEAM color scale helper
- [`scale_gleam_fill()`](https://JamesWu7.github.io/GLEAM/reference/scale_gleam_fill.md)
  : GLEAM fill scale helper

## Export and method comparison

- [`collect_scores()`](https://JamesWu7.github.io/GLEAM/reference/collect_scores.md)
  : Collect score objects from multiple methods
- [`summarize_scores()`](https://JamesWu7.github.io/GLEAM/reference/summarize_scores.md)
  : Summarize score table
- [`pivot_scores_long()`](https://JamesWu7.github.io/GLEAM/reference/pivot_scores_long.md)
  : Pivot score matrix to long format
- [`export_scores()`](https://JamesWu7.github.io/GLEAM/reference/export_scores.md)
  : Export score table

## Example data

- [`hallmark`](https://JamesWu7.github.io/GLEAM/reference/hallmark.md) :
  Built-in hallmark-like gene sets
- [`immune_small`](https://JamesWu7.github.io/GLEAM/reference/immune_small.md)
  : Built-in immune-focused small gene sets
- [`pbmc_medium_matrix`](https://JamesWu7.github.io/GLEAM/reference/pbmc_medium_matrix.md)
  : Medium PBMC expression matrix for realistic demos
- [`pbmc_medium_meta`](https://JamesWu7.github.io/GLEAM/reference/pbmc_medium_meta.md)
  : Medium PBMC metadata for realistic demos
- [`pbmc_small`](https://JamesWu7.github.io/GLEAM/reference/pbmc_small.md)
  : Small PBMC example dataset (legacy compatibility)
- [`spatial_medium_coords`](https://JamesWu7.github.io/GLEAM/reference/spatial_medium_coords.md)
  : Medium spatial coordinates
- [`spatial_medium_expr`](https://JamesWu7.github.io/GLEAM/reference/spatial_medium_expr.md)
  : Medium spatial expression matrix
- [`spatial_medium_meta`](https://JamesWu7.github.io/GLEAM/reference/spatial_medium_meta.md)
  : Medium spatial metadata
- [`toy_expr`](https://JamesWu7.github.io/GLEAM/reference/toy_expr.md) :
  Tiny toy expression dataset for tests
