# GLEAM API Migration: pathway -> signature (breaking cleanup)

GLEAM now uses **signature** as the canonical umbrella term in the public API.
Legacy `*_pathway*` names are no longer exported in v0.2+.
Internal compatibility wrappers remain in code for controlled migration, but are removed from public reference pages.

## Canonical names

- `score_signature()`
- `aggregate_signature()`
- `test_signature()`
- `test_signature_spatial()`
- `test_signature_trajectory()`

## Workflow entry points

- `run_gleam()`

## Public API classification (v0.2.0)

- Core analysis API:
  - `score_signature()`
  - `aggregate_signature()`
  - `test_signature()`
  - `run_gleam()`
- Scenario-specific testing:
  - `test_signature_spatial()`
  - `test_signature_trajectory()`
- Trajectory/spatial utilities:
  - `as_trajectory_data()`, `extract_pseudotime()`, `extract_lineage()`
  - `join_score_trajectory()`, `map_scores_to_trajectory()`, `join_score_spatial()`
- Visualization API (kept separate by design):
  - `plot_embedding_score()`, `plot_spatial_score()`, `plot_spatial_multi()`, `plot_spatial_compare()`
  - `plot_dot()`, `plot_dot_bar()`, `plot_violin()`, `plot_split_violin()`, `plot_ridge()`
  - `plot_box()`, `plot_heatmap()`, `plot_volcano()`, `plot_pseudotime_score()`, `plot_trajectory_score()`
  - `plot_celltype_compare()`, `plot_group_in_celltype()`, `plot_pseudobulk_box()`
- Geneset/signature management:
  - `get_geneset()`, `list_geneset_sources()`, `search_geneset()`
  - `as_geneset()`, `read_gmt()`, `check_geneset()`, `match_geneset()`
- Export/summary/comparison:
  - `collect_scores()`, `summarize_scores()`, `pivot_scores_long()`, `export_scores()`, `compare_scoring_methods()`
- Theme/palette helpers:
  - `gleam_theme()`, `apply_gleam_theme()`, `get_palette()`, `list_palettes()`, `scale_gleam_color()`, `scale_gleam_fill()`

## Notes

- Plotting functions remain separate by design.
- Differential entrypoints were merged into the `test_signature*` family.
- Use `compare_scoring_methods()` (the `compare_methods()` alias is removed from exports).
