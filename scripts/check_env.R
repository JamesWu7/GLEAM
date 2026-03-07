message("[GLEAM] Environment check")
message("R version: ", R.version.string)

runtime <- c("Matrix", "ggplot2")
dev <- c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown", "pkgload", "usethis", "remotes")
optional <- c(
  "Seurat", "SeuratObject", "msigdbr", "UCell", "AUCell", "GSVA", "singscore",
  "slingshot", "monocle3", "tradeSeq", "patchwork", "viridisLite", "RColorBrewer"
)

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
check_pkgs(optional, "Optional integration packages")

if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("Seurat", quietly = TRUE)) {
  message("\n[GLEAM] Note: Seurat mode (seurat = TRUE) requires Seurat/SeuratObject.")
}

message("\n[GLEAM] Done.")
