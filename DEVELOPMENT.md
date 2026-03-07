# GLEAM Development Guide

## Fresh environment bootstrap
```r
source("scripts/setup_dev.R")
source("scripts/check_env.R")
source("scripts/check_optional_deps.R")
```

## Core development commands
```r
devtools::document()
devtools::test()
devtools::check()
devtools::build()
```

## Refresh example data
```r
source("scripts/refresh_example_data.R")
```

## Build example Seurat object
```r
source("scripts/build_example_seurat.R")
```

## Render tutorials to HTML
```r
source("scripts/render_examples.R")
```

Rendered files are placed in `docs/tutorial_html/`.

## Notes
- Matrix workflow is first-class and does not require Seurat.
- Seurat workflows support v4 slot and v5 layer conventions.
- Optional integrations are dependency-gated with explicit install messages.
