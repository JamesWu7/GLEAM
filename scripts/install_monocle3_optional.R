cat("[GLEAM] Optional monocle3 installation helper\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("[GLEAM] Installing BiocManager...\n")
  install.packages("BiocManager")
}

cat("[GLEAM] Installing monocle3 from Bioconductor...\n")
BiocManager::install("monocle3", ask = FALSE, update = FALSE)

ok <- requireNamespace("monocle3", quietly = TRUE)
cat("[GLEAM] monocle3 status: ", ifelse(ok, "OK", "NOT INSTALLED"), "\n", sep = "")
if (!ok) {
  cat("[GLEAM] monocle3 remains optional. Core GLEAM workflows do not require it.\n")
}
