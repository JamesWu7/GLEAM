message("[scPathway] Environment check")
message("R version: ", R.version.string)

runtime <- c("Matrix", "ggplot2")
dev <- c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "pkgload", "usethis", "remotes")
optional <- c("Seurat", "SeuratObject")

check_pkgs <- function(pkgs, label) {
  message("\n", label)
  for (pkg in pkgs) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    status <- if (ok) "OK" else "MISSING"
    message(sprintf("  - %-12s : %s", pkg, status))
  }
}

check_pkgs(runtime, "Runtime packages")
check_pkgs(dev, "Development packages")
check_pkgs(optional, "Optional Seurat compatibility packages")

if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
  message("\n[scPathway] Note: Seurat mode (seurat = TRUE) will require Seurat or SeuratObject.")
}

message("\n[scPathway] Done.")
