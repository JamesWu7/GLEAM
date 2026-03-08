# GLEAM Example Data Layout

## Data tiers

1. Packaged lightweight data (`data/*.rda`)
- Used for routine examples and unit tests.
- Includes matrix/meta toy/medium datasets and genesets.

2. Full Seurat demo objects (`inst/extdata/full_examples/*.rds`)
- Used for realistic documentation workflows.
- Excluded from package build via `.Rbuildignore` to avoid package bloat.

3. Lightweight Seurat test subsets (`inst/extdata/test_examples/*.rds`)
- Derived from full objects using `data-raw/derive_test_subsets.R`.
- Intended for optional integration tests and CI smoke checks.
