cat("[GLEAM] Optional monocle3 installation helper\n")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

cat("[GLEAM] Ensuring BiocManager is available...\n")
install_if_missing("BiocManager")

bioc_version <- Sys.getenv("GLEAM_BIOC_VERSION", unset = "")
if (nzchar(bioc_version)) {
  cat("[GLEAM] Setting Bioconductor version to ", bioc_version, " (optional override)\n", sep = "")
  BiocManager::install(version = bioc_version, ask = FALSE, update = FALSE)
} else {
  cat("[GLEAM] Using default Bioconductor version for your R release.\n")
}

bioc_pkgs <- c(
  "BiocGenerics", "DelayedArray", "DelayedMatrixStats",
  "limma", "lme4", "S4Vectors", "SingleCellExperiment",
  "SummarizedExperiment", "batchelor", "HDF5Array",
  "ggrastr"
)
cat("[GLEAM] Installing key Bioconductor dependencies for monocle3...\n")
BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

cat("[GLEAM] Installing remotes for GitHub package install...\n")
install_if_missing("remotes")

cat("[GLEAM] Installing BPCells from GitHub (optional monocle3 ecosystem dependency)...\n")
remotes::install_github("bnprks/BPCells/r", upgrade = "never", dependencies = FALSE)

cat("[GLEAM] Installing monocle3 from Bioconductor...\n")
BiocManager::install("monocle3", ask = FALSE, update = FALSE)

ok <- requireNamespace("monocle3", quietly = TRUE)
cat("[GLEAM] monocle3 status: ", ifelse(ok, "OK", "NOT INSTALLED"), "\n", sep = "")
if (!ok) {
  cat("[GLEAM] monocle3 remains optional. Core GLEAM workflows do not require it.\n")
  cat("[GLEAM] Tip: you can rerun this script with GLEAM_BIOC_VERSION=3.21 if needed.\n")
}
