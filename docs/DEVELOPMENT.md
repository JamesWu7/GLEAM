# GLEAM Development Guide

## Fresh environment bootstrap

``` r
source("scripts/setup_dev.R")
source("scripts/check_env.R")
source("scripts/check_optional_deps.R")
source("scripts/check_monocle3.R")
```

## Core development commands

``` r
devtools::document()
devtools::test()
devtools::check()
devtools::build()
```

## Refresh example data

``` r
source("scripts/refresh_example_data.R")
```

## Build example Seurat object

``` r
source("scripts/build_example_seurat.R")
```

## Generate homepage workflow figures

``` r
source("scripts/generate_homepage_figures.R")
```

Generated figures are placed in `assets/figures/` and copied into
`docs/assets/` during pkgdown deployment.

## Audit full Seurat example object structure (optional)

``` r
source("scripts/audit_example_objects.R")
```

## Optional Monocle3 helper

``` r
source("scripts/install_monocle3_optional.R")
```

## Build pkgdown site

``` r
pkgdown::build_site()
```

## Citation workflow

``` r
citation("GLEAM")
```

Citation metadata is maintained in `inst/CITATION`.

## Notes

- Matrix workflow is first-class and does not require Seurat.
- Seurat workflows support v4 slot and v5 layer conventions.
- Optional integrations are dependency-gated with explicit install
  messages.
