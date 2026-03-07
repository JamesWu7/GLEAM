# Contributing to scPathway

## Development setup
1. `source("scripts/setup_dev.R")`
2. `source("scripts/check_env.R")`
3. `devtools::load_all()`

## Development cycle
1. Create a feature branch from `main`.
2. Add or update roxygen documentation.
3. Add or update tests under `tests/testthat`.
4. Run:
   - `devtools::document()`
   - `devtools::test()`
   - `devtools::check()`
5. Update `NEWS.md` for user-visible changes.

## Coding standards
- Keep runtime dependencies minimal.
- Prefer matrix-native internals.
- Keep Seurat support optional.
- Add explicit error and warning messages.
