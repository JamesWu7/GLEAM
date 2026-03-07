pkgs <- c(
  Seurat = "Seurat v4/v5 input",
  SeuratObject = "Seurat object access",
  msigdbr = "MSigDB/GO/KEGG/Reactome geneset sources",
  UCell = "UCell scoring wrapper",
  AUCell = "AUCell scoring wrapper",
  GSVA = "GSVA scoring wrapper",
  singscore = "singscore wrapper",
  slingshot = "slingshot trajectory interface",
  monocle3 = "monocle3 trajectory interface",
  tradeSeq = "trajectory differential testing",
  ggridges = "ridge visualization",
  patchwork = "multi-panel composition"
)

cat("Optional dependency report\n")
for (p in names(pkgs)) {
  ok <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("- %-12s : %-7s (%s)\n", p, if (ok) "OK" else "MISSING", pkgs[[p]]))
}
