# GLEAM API Migration: pathway -\> signature (breaking cleanup)

GLEAM now uses **signature** as the canonical umbrella term in the
public API. Legacy `*_pathway*` names are no longer exported in v0.2+.
Internal compatibility wrappers remain in code for controlled migration,
but are removed from public reference pages.

## Canonical names

- [`score_signature()`](https://JamesWu7.github.io/GLEAM/reference/score_signature.md)
- [`aggregate_signature()`](https://JamesWu7.github.io/GLEAM/reference/aggregate_signature.md)
- [`test_signature()`](https://JamesWu7.github.io/GLEAM/reference/test_signature.md)
- [`test_signature_spatial()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_spatial.md)
- [`test_signature_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_trajectory.md)

## Workflow entry points

- [`run_gleam()`](https://JamesWu7.github.io/GLEAM/reference/run_gleam.md)

## Public API classification (v0.2.0)

- Core analysis API:
  - [`score_signature()`](https://JamesWu7.github.io/GLEAM/reference/score_signature.md)
  - [`aggregate_signature()`](https://JamesWu7.github.io/GLEAM/reference/aggregate_signature.md)
  - [`test_signature()`](https://JamesWu7.github.io/GLEAM/reference/test_signature.md)
  - [`run_gleam()`](https://JamesWu7.github.io/GLEAM/reference/run_gleam.md)
- Scenario-specific testing:
  - [`test_signature_spatial()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_spatial.md)
  - [`test_signature_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/test_signature_trajectory.md)
- Trajectory/spatial utilities:
  - [`as_trajectory_data()`](https://JamesWu7.github.io/GLEAM/reference/as_trajectory_data.md),
    [`extract_pseudotime()`](https://JamesWu7.github.io/GLEAM/reference/extract_pseudotime.md),
    [`extract_lineage()`](https://JamesWu7.github.io/GLEAM/reference/extract_lineage.md)
  - [`join_score_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/join_score_trajectory.md),
    [`map_scores_to_trajectory()`](https://JamesWu7.github.io/GLEAM/reference/map_scores_to_trajectory.md),
    [`join_score_spatial()`](https://JamesWu7.github.io/GLEAM/reference/join_score_spatial.md)
- Visualization API (kept separate by design):
  - [`plot_embedding_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_embedding_score.md),
    [`plot_spatial_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_score.md),
    [`plot_spatial_multi()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_multi.md),
    [`plot_spatial_compare()`](https://JamesWu7.github.io/GLEAM/reference/plot_spatial_compare.md)
  - [`plot_dot()`](https://JamesWu7.github.io/GLEAM/reference/plot_dot.md),
    [`plot_dot_bar()`](https://JamesWu7.github.io/GLEAM/reference/plot_dot_bar.md),
    [`plot_violin()`](https://JamesWu7.github.io/GLEAM/reference/plot_violin.md),
    [`plot_split_violin()`](https://JamesWu7.github.io/GLEAM/reference/plot_split_violin.md),
    [`plot_ridge()`](https://JamesWu7.github.io/GLEAM/reference/plot_ridge.md)
  - [`plot_box()`](https://JamesWu7.github.io/GLEAM/reference/plot_box.md),
    [`plot_heatmap()`](https://JamesWu7.github.io/GLEAM/reference/plot_heatmap.md),
    [`plot_volcano()`](https://JamesWu7.github.io/GLEAM/reference/plot_volcano.md),
    [`plot_pseudotime_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_pseudotime_score.md),
    [`plot_trajectory_score()`](https://JamesWu7.github.io/GLEAM/reference/plot_trajectory_score.md)
  - [`plot_celltype_compare()`](https://JamesWu7.github.io/GLEAM/reference/plot_celltype_compare.md),
    [`plot_group_in_celltype()`](https://JamesWu7.github.io/GLEAM/reference/plot_group_in_celltype.md),
    [`plot_pseudobulk_box()`](https://JamesWu7.github.io/GLEAM/reference/plot_pseudobulk_box.md)
- Geneset/signature management:
  - [`get_geneset()`](https://JamesWu7.github.io/GLEAM/reference/get_geneset.md),
    [`list_geneset_sources()`](https://JamesWu7.github.io/GLEAM/reference/list_geneset_sources.md),
    [`search_geneset()`](https://JamesWu7.github.io/GLEAM/reference/search_geneset.md)
  - [`as_geneset()`](https://JamesWu7.github.io/GLEAM/reference/as_geneset.md),
    [`read_gmt()`](https://JamesWu7.github.io/GLEAM/reference/read_gmt.md),
    [`check_geneset()`](https://JamesWu7.github.io/GLEAM/reference/check_geneset.md),
    [`match_geneset()`](https://JamesWu7.github.io/GLEAM/reference/match_geneset.md)
- Export/summary/comparison:
  - [`collect_scores()`](https://JamesWu7.github.io/GLEAM/reference/collect_scores.md),
    [`summarize_scores()`](https://JamesWu7.github.io/GLEAM/reference/summarize_scores.md),
    [`pivot_scores_long()`](https://JamesWu7.github.io/GLEAM/reference/pivot_scores_long.md),
    [`export_scores()`](https://JamesWu7.github.io/GLEAM/reference/export_scores.md),
    [`compare_scoring_methods()`](https://JamesWu7.github.io/GLEAM/reference/compare_scoring_methods.md)
- Theme/palette helpers:
  - [`gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/gleam_theme.md),
    [`apply_gleam_theme()`](https://JamesWu7.github.io/GLEAM/reference/apply_gleam_theme.md),
    [`get_palette()`](https://JamesWu7.github.io/GLEAM/reference/get_palette.md),
    [`list_palettes()`](https://JamesWu7.github.io/GLEAM/reference/list_palettes.md),
    [`scale_gleam_color()`](https://JamesWu7.github.io/GLEAM/reference/scale_gleam_color.md),
    [`scale_gleam_fill()`](https://JamesWu7.github.io/GLEAM/reference/scale_gleam_fill.md)

## Notes

- Plotting functions remain separate by design.
- Differential entrypoints were merged into the `test_signature*`
  family.
- Use
  [`compare_scoring_methods()`](https://JamesWu7.github.io/GLEAM/reference/compare_scoring_methods.md)
  (the `compare_methods()` alias is removed from exports).
