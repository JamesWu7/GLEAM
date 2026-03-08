# GLEAM Repository Diagnosis and Implementation Plan (Phase A Baseline)

Date: 2026-03-08
Repository: `JamesWu7/GLEAM`

## 1) Gap analysis

### Already present
- Matrix-first scoring/testing workflow and score/test object classes (`gleam_score`, `gleam_test`).
- Seurat input path (v4/v5 oriented), spatial utilities, and trajectory utility entry points.
- Differential APIs (`test_*`, `differential_*`, `compare_*`) and plotting family (`plot_*`) already separated by task.
- Multiple vignettes and pkgdown scaffolding present.

### Main gaps / breakpoints found
- Canonical naming was still pathway-first across public API and docs.
- Monocle3 optional backend boundary was not explicit enough (dependency checks dispersed).
- Missing centralized optional dependency assertions and user-facing install hints.
- Stale repository owner links (`jameswoo`) in metadata/docs.
- README/docs presentation not sufficiently figure-first; docs URL currently not live (depends on workflow deploy).
- Full demo Seurat objects were root-scattered before refactor.
- CI lacked explicit split between core and optional Monocle3 backend validation.
- Exported helper surface had mixed abstraction levels and unclear canonical-vs-compatibility positioning.

## 2) Current trajectory architecture summary
- Public trajectory entry points: `as_trajectory_data()`, `test_pathway_trajectory()`, `plot_pseudotime_score()`, `plot_trajectory_score()`, `trajectory_table()`.
- Inputs can be vectors/metadata/object carriers; extraction existed but backend-specific handling was interleaved.
- Current refactor introduces explicit backend path (`backend = c("auto", "monocle3", "internal")`) and normalized result conversion.

## 3) Dependency risk analysis
- Optional ecosystem packages are numerous; risk is accidental hard dependency in code paths/examples/tests.
- Highest risk package: `monocle3` (non-core trajectory backend with fragile install paths).
- Seurat-related risks are contained by `requireNamespace` checks and test skips.
- pkgdown local build currently blocked in this environment by `systemfonts/textshaping/ragg` toolchain.

## 4) DESCRIPTION / Suggests / Imports / Remotes audit
- `monocle3` retained in `Suggests` (not in `Imports`/`Depends`) — correct.
- Added `Biobase` to `Suggests` to match `Biobase::pData` usage in trajectory utilities.
- `Matrix` remains in `Imports`; now used explicitly via `Matrix::colMeans` in scoring helper path.
- URL/BugReports corrected to `https://github.com/JamesWu7/GLEAM`.

## 5) Test structure audit
- Core tests now run without Monocle3.
- Optional Monocle3 tests introduced and gated by `skip_if_not_installed("monocle3")`.
- Added compatibility tests for signature aliases and legacy pathway wrappers.

## 6) Vignette / tutorial / README audit
- Main vignettes are `.Rmd` and now include chunk sizing controls.
- Seurat tutorials show object context early (`dim`, object print, metadata preview, UMAP).
- README now uses figure-first sections and signature-first API examples.
- Removed accidental generated vignette HTML artifact from `vignettes/`.

## 7) CI risk audit
- Core workflow: `.github/workflows/R-CMD-check.yaml` with `_R_CHECK_FORCE_SUGGESTS_=false`.
- Optional workflow added: `.github/workflows/optional-monocle3.yaml`.
- pkgdown workflow builds figures then site.
- Remaining risk: monocle3 install reliability on GitHub runner matrix.

## 8) Public API analysis
- Canonical family now introduced: `score_signature()`, `aggregate_signature()`, `test_signature()` and spatial/trajectory variants.
- Legacy pathway family retained as compatibility wrappers.
- `run_gleam()` canonical; `run_scpathway()` retained with deprecation warning.

## 9) Result-object analysis
- Score object: stable (`score`, `meta`, `method`, geneset metadata, params).
- Test object: stable table-centric output with comparison metadata.
- New trajectory result standardization added via `new_gleam_trajectory_result()` and `as_trajectory_result()`.

## 10) Public function inventory (exported)

### Core workflow
- `run_gleam`, `run_scpathway`
- `score_signature`, `aggregate_signature`, `test_signature`, `differential_signature`
- `score_pathway`, `aggregate_pathway`, `test_pathway`, `differential_pathway`
- `compare_celltypes`, `compare_groups_within_celltype`

### Geneset
- `get_geneset`, `list_geneset_sources`, `search_geneset`, `as_geneset`, `read_gmt`, `check_geneset`, `match_geneset`

### Trajectory
- `as_trajectory_data`, `extract_pseudotime`, `extract_lineage`, `join_score_trajectory`, `map_scores_to_trajectory`, `trajectory_table`
- `test_signature_trajectory`, `differential_signature_trajectory`
- `test_pathway_trajectory`, `differential_pathway_trajectory`
- `plot_pseudotime_score`, `plot_trajectory_score`

### Spatial
- `join_score_spatial`, `spatial_table`
- `test_signature_spatial`, `differential_signature_spatial`
- `test_pathway_spatial`, `differential_pathway_spatial`
- `plot_spatial_score`, `plot_spatial_compare`, `plot_spatial_multi`

### Visualization
- `plot_violin`, `plot_split_violin`, `plot_ridge`, `plot_box`, `plot_dot`, `plot_dot_bar`, `plot_heatmap`, `plot_volcano`, `plot_embedding_score`, `plot_celltype_compare`, `plot_group_in_celltype`, `plot_pseudobulk_box`
- `gleam_theme`, `apply_gleam_theme`
- `get_palette`, `list_palettes`, `scale_gleam_color`, `scale_gleam_fill`

### Export/comparison
- `collect_scores`, `summarize_scores`, `pivot_scores_long`, `export_scores`, `compare_scoring_methods`, `compare_methods`

### Other helpers currently exported
- `extract_embedding`, `seurat_mode`

## 11) Naming consistency audit
- Resolved: canonical shift to signature family for major analysis APIs.
- Retained wrappers: pathway family kept for backward compatibility.
- Remaining ambiguity (intended): `test_*` and `differential_*` both exist; differential is alias-oriented for user familiarity.
- Plot APIs intentionally remain separate (no monolithic plot API merge).
- `seurat_mode` remains public for now; candidate for future demotion.

### Explicit pathway -> signature mapping

| Legacy name | Canonical name | Decision | Compatibility plan |
|---|---|---|---|
| `score_pathway` | `score_signature` | Rename and keep public | keep `score_pathway()` wrapper and docs alias |
| `aggregate_pathway` | `aggregate_signature` | Rename and keep public | keep wrapper + docs alias |
| `test_pathway` | `test_signature` | Rename and keep public | keep wrapper + docs alias |
| `differential_pathway` | `differential_signature` | Rename and keep public | keep wrapper; internally alias to test family |
| `test_pathway_spatial` | `test_signature_spatial` | Rename and keep public | keep wrapper + docs alias |
| `differential_pathway_spatial` | `differential_signature_spatial` | Rename and keep public | keep wrapper alias |
| `test_pathway_trajectory` | `test_signature_trajectory` | Rename and keep public | keep wrapper + backend-aware forwarding |
| `differential_pathway_trajectory` | `differential_signature_trajectory` | Rename and keep public | keep wrapper alias |
| `run_scpathway` | `run_gleam` | Deprecate old brand | keep wrapper with warning, remove from canonical examples |

## 12) Pathway vs signature terminology audit
- Canonical docs/examples now prioritize `score_signature()/test_signature()`.
- Pathway terms preserved where semantically specific (external pathway collections, compatibility wrappers, table columns).

## 13) Vignette rendering / figure reuse strategy audit
- Vignette Rmds tuned for readable HTML output sizing.
- Stable homepage/doc figure export script: `scripts/generate_homepage_figures.R`.
- Figure targets in `docs/assets/figures/` and surfaced in README.

## 14) README/docs landing visual audit
- README top now has larger logo and figure-first snapshot block.
- Navbar-style links for Documentation/Reference/Tutorials/Citation added.
- Remaining deployment dependency: GitHub Pages publish run must succeed for live docs URL.

## 15) Example data / asset placement audit
- Logo moved to stable docs path: `docs/assets/GLEAM_LOG.jpg`.
- Full example Seurat objects moved from repo root to `inst/extdata/full_examples/`.
- Packaged lightweight matrix/meta data retained in `data/`.

## 16) Example data size / test-subset strategy audit
- Full objects retained for realistic docs/demo.
- Test-speed matrix datasets in `data/` remain primary for routine tests.
- Subset derivation script added: `data-raw/derive_test_subsets.R` (outputs to `inst/extdata/test_examples/` when Seurat is available).

## 17) Proposed backend boundary design (trajectory)
- Input conversion: `.as_monocle3_input()`.
- Backend execution: `.run_trajectory_monocle3()`.
- Result normalization: `new_gleam_trajectory_result()`, `as_trajectory_result()`.

## 18) Proposed dependency-gating design
- Centralized optional checks in `R/optional_deps.R`.
- Dedicated `.assert_monocle3()` with concise optional-install guidance.
- Route all optional backend paths through shared assertions.

## 19) Proposed test split strategy
- Core tests: default `testthat` suite independent of Monocle3.
- Optional tests: explicit `trajectory-monocle3` / optional dependency behavior tests.

## 20) Proposed CI split strategy
- Core: `R-CMD-check.yaml` (no Monocle3 requirement).
- Optional: `optional-monocle3.yaml` installs monocle3 and runs backend-specific tests.
- Website: `pkgdown.yaml` generates figures + site.

## 21) Proposed documentation strategy
- Signature-first examples in README and vignettes.
- Explicit install tiers (core/recommended/optional trajectory).
- Dedicated migration/integration docs: `API_MIGRATION.md`, `MONOCLE3_INTEGRATION.md`.

## 22) Install guidance strategy
- Keep package install lightweight.
- Optional backend install documented with explicit user guidance and scripts.
- No auto-install behavior in package load path.

## 23) API cleanup strategy
- Canonical entry points: signature family + `run_gleam`.
- Compatibility wrappers retained with migration docs.
- Keep plotting functions separate and discoverable.

## 24) Deprecation / compatibility strategy
- `run_scpathway()` kept as deprecating wrapper to `run_gleam()`.
- Pathway family kept as wrappers around signature-preferred docs.
- `compare_methods()` retained as alias to `compare_scoring_methods()`.

## 25) Vignette knitting + figure surfacing strategy
- Vignette chunks include explicit sizing (`fig.width`, `fig.height`, `out.width`, `dpi`).
- Reusable figures generated from package workflows into `docs/assets/figures`.
- README and pkgdown home consume stable figure paths.

## 26) Example-data and asset layout strategy
- `data/`: packaged lightweight test/demo matrices/meta.
- `inst/extdata/full_examples/`: larger full Seurat objects.
- `data-raw/`: derivation scripts and data process notes.
- `docs/assets/`: logo + presentation figures.

## 27) Subset-for-testing strategy
- Keep core tests on lightweight packaged matrix data.
- Optionally produce Seurat subset RDS files for CI developer workflows using `data-raw/derive_test_subsets.R`.

## 28) Homepage visual strategy (irGSEA/liana-inspired principles)
- Strong top hierarchy: logo -> title/subtitle -> nav links -> workflow figures.
- Figure-first evidence of scRNA/spatial/trajectory outputs near top.
- Avoid text-only scaffold presentation.

## 29) File-by-file implementation plan (current + next)
- Dependency hardening: `R/optional_deps.R`, `R/utils_general.R`, `DESCRIPTION`.
- Signature API: `R/api_signature.R`, docs aliases in `man/*pathway*.Rd`, `NAMESPACE`.
- Trajectory backend: `R/trajectory_backend_monocle3.R`, `R/trajectory_result.R`, `R/trajectory_utils.R`, `R/api_test.R`, `R/api_run.R`, `R/input_trajectory.R`.
- Tests: new trajectory optional/compatibility tests under `tests/testthat/`.
- CI: `.github/workflows/optional-monocle3.yaml`, updated `pkgdown.yaml`.
- Docs/homepage: `README.md`, `README.Rmd`, `_pkgdown.yml`, vignettes, figure script.
- Data/assets: `inst/extdata/full_examples/`, `data-raw/derive_test_subsets.R`, `docs/assets/GLEAM_LOG.jpg`.

## 30) Migration / backward-compatibility plan
- New users: use signature-first API.
- Existing users: pathway wrappers continue to work.
- Migration details documented in `API_MIGRATION.md`.

## 31) Repository metadata consistency plan
- All owner/repo links standardized to `JamesWu7/GLEAM`.
- pkgdown URL standardized to `https://jameswu7.github.io/GLEAM`.

## 32) Work-tree / task-tree for modular execution

### Module tree
1. Metadata + dependency hardening
2. Signature API migration + compatibility wrappers
3. Trajectory optional backend boundary
4. Spatial integration stability
5. Visualization/readability layer
6. Vignette and figure generation
7. Docs publishing + homepage UX
8. CI split and test gating

### Task dependency order
1 -> 2 -> 3 -> 8 -> 6 -> 7
4 and 5 run in parallel after 2; feed into 6 and 7.

### Code/docs/data/test/publish split
- Code: `R/`
- Tests: `tests/testthat/`
- Data: `data/`, `inst/extdata/`, `data-raw/`
- Documentation: `README*`, `vignettes/`, `man/`, `API_MIGRATION.md`, `MONOCLE3_INTEGRATION.md`
- Publishing: `_pkgdown.yml`, `.github/workflows/pkgdown.yaml`, `docs/assets/`

## 33) Exported function classification table

### Keep as core public API
- `score_signature`, `aggregate_signature`, `test_signature`, `differential_signature`
- `run_gleam`
- `get_geneset`, `list_geneset_sources`, `search_geneset`, `as_geneset`, `read_gmt`, `check_geneset`, `match_geneset`
- `test_signature_spatial`, `test_signature_trajectory`
- plotting family: `plot_violin`, `plot_split_violin`, `plot_ridge`, `plot_box`, `plot_dot`, `plot_dot_bar`, `plot_heatmap`, `plot_volcano`, `plot_embedding_score`, `plot_celltype_compare`, `plot_group_in_celltype`, `plot_pseudobulk_box`, `plot_spatial_score`, `plot_spatial_compare`, `plot_spatial_multi`, `plot_pseudotime_score`, `plot_trajectory_score`

### Convenience wrappers (keep public)
- `compare_celltypes`, `compare_groups_within_celltype`
- `differential_signature_spatial`, `differential_signature_trajectory`
- `compare_methods` (alias)

### Lower-level utilities (keep public)
- `as_trajectory_data`, `extract_pseudotime`, `extract_lineage`, `join_score_trajectory`, `map_scores_to_trajectory`, `trajectory_table`
- `join_score_spatial`, `spatial_table`
- `collect_scores`, `summarize_scores`, `pivot_scores_long`, `export_scores`, `compare_scoring_methods`
- palette/theme helpers and `extract_embedding`

### Deprecated compatibility wrappers (keep for migration)
- `score_pathway`, `aggregate_pathway`, `test_pathway`, `differential_pathway`
- `test_pathway_spatial`, `differential_pathway_spatial`
- `test_pathway_trajectory`, `differential_pathway_trajectory`
- `run_scpathway`

### Candidate to internalize in next major release
- `seurat_mode` (currently exported for backward compatibility; role is mostly internal mode inference)

## 34) Risk register
- R1: Monocle3 install instability on CI runners.
  - Mitigation: separate optional CI + clear skips.
- R2: Docs URL remains 404 until gh-pages deployment succeeds.
  - Mitigation: pkgdown workflow with deploy step + debug after first run.
- R3: Local pkgdown build blocked by system font toolchain.
  - Mitigation: rely on CI build environment; keep figure generation script independent.
- R4: API confusion from dual pathway/signature naming.
  - Mitigation: signature-first docs + migration page + compatibility tests.
- R5: Large full example data inflating package footprint if accidentally included.
  - Mitigation: `.Rbuildignore` excludes `inst/extdata/full_examples`.
