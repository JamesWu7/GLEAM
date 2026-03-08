# Monocle3 Integration in GLEAM

Monocle3 support in GLEAM is optional.

## Design

1. Core GLEAM scoring/testing workflows do not require Monocle3.
2. Monocle3 is declared in `Suggests`, not `Imports`/`Depends`.
3. Monocle3 code paths are runtime-gated through optional dependency checks.
4. Trajectory utilities normalize backend outputs into a GLEAM trajectory result structure.

## When Monocle3 is required

- Only for Monocle3-specific trajectory backend usage (for example `backend = "monocle3"` in trajectory testing).
- Not required for matrix scoring, differential analysis, spatial testing, or plotting in non-trajectory mode.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Optional: pin Bioconductor version when troubleshooting
# BiocManager::install(version = "3.21")

BiocManager::install(c(
  "BiocGenerics", "DelayedArray", "DelayedMatrixStats",
  "limma", "lme4", "S4Vectors", "SingleCellExperiment",
  "SummarizedExperiment", "batchelor", "HDF5Array",
  "ggrastr"
))

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("bnprks/BPCells/r")

BiocManager::install("monocle3")
```

Notes:

- `BiocManager::install(version = "3.21")` is a troubleshooting option, not a mandatory fixed version.
- If installation fails, you can continue using non-trajectory GLEAM workflows.
- The maintainer helper script `scripts/install_monocle3_optional.R` automates this optional path.
