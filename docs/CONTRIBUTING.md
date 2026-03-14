# Contributing to GLEAM

## Development setup

1.  `source("scripts/setup_dev.R")`
2.  `source("scripts/check_env.R")`
3.  [`devtools::load_all()`](https://devtools.r-lib.org/reference/load_all.html)

## Development cycle

1.  Create a feature branch from `main`.
2.  Update roxygen documentation.
3.  Add/update tests under `tests/testthat`.
4.  Run:
    - [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
    - [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
    - [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
5.  Update `NEWS.md` for user-visible changes.

## Principles

- Keep runtime dependencies minimal.
- Keep matrix workflows robust.
- Keep Seurat/trajectory/spatial integrations optional and explicit.
- Use clear warnings and errors for optional feature gates.
