# Development Guide

## Fresh environment bootstrap
```r
source("scripts/setup_dev.R")
source("scripts/check_env.R")
```

## Core development commands
```r
devtools::document()
devtools::test()
devtools::check()
devtools::build()
```

## Notes
- Seurat v5 support is optional and only required when `seurat = TRUE`.
- Matrix workflow is fully supported without Seurat.
