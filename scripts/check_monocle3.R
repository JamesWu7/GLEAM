cat("[GLEAM] Optional Monocle3 environment check\n")
cat("R version: ", R.version.string, "\n", sep = "")

bioc_ver <- tryCatch({
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    as.character(BiocManager::version())
  } else {
    NA_character_
  }
}, error = function(e) NA_character_)
cat("Bioconductor version: ", ifelse(is.na(bioc_ver), "not detected", bioc_ver), "\n", sep = "")

has_bioc <- requireNamespace("BiocManager", quietly = TRUE)
has_m3 <- requireNamespace("monocle3", quietly = TRUE)

cat("BiocManager: ", ifelse(has_bioc, "OK", "MISSING"), "\n", sep = "")
cat("monocle3: ", ifelse(has_m3, "OK", "MISSING"), "\n", sep = "")

if (!has_bioc) {
  cat("Hint: install.packages('BiocManager')\n")
}
if (!has_m3) {
  cat("Hint: BiocManager::install('monocle3')\n")
  cat("Note: monocle3 is optional and only needed for trajectory workflows.\n")
}
