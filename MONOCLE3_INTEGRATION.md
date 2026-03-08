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
install.packages("BiocManager")
BiocManager::install("monocle3")
```

If installation fails, you can continue using non-trajectory GLEAM workflows.
