message("[scPathway] Checking development dependencies...")

required <- c(
  "devtools", "roxygen2", "testthat", "knitr",
  "rmarkdown", "pkgload", "usethis", "remotes"
)

is_installed <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

missing <- required[!vapply(required, is_installed, logical(1))]

if (length(missing) == 0) {
  message("[scPathway] All development dependencies are already installed.")
} else {
  message("[scPathway] Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE)
  message("[scPathway] Installation complete.")
}

message("[scPathway] Next steps:")
message("  source('scripts/check_env.R')")
message("  devtools::document()")
message("  devtools::test()")
message("  devtools::check()")
