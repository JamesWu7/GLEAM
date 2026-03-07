message("[GLEAM] Checking development dependencies...")

required <- c(
  "devtools", "roxygen2", "testthat", "knitr",
  "rmarkdown", "pkgload", "usethis", "remotes"
)

optional <- c(
  "Seurat", "SeuratObject", "msigdbr", "patchwork", "viridisLite"
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

missing <- required[!vapply(required, is_installed, logical(1))]

if (length(missing) == 0) {
  message("[GLEAM] All development dependencies are already installed.")
} else {
  message("[GLEAM] Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
  message("[GLEAM] Installation complete.")
}

missing_optional <- optional[!vapply(optional, is_installed, logical(1))]
if (length(missing_optional) > 0) {
  message("[GLEAM] Optional packages not installed: ", paste(missing_optional, collapse = ", "))
}

message("[GLEAM] Next steps:")
message("  source('scripts/check_env.R')")
message("  devtools::document()")
message("  devtools::test()")
message("  devtools::check()")
message("  source('scripts/render_examples.R')")
