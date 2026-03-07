set.seed(1001)

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required to generate package data.")
}

# Built-in genesets
immune_small <- list(
  T_cell_activation = c("CD3D", "CD3E", "LCK", "TRAC", "IL7R"),
  Cytotoxicity = c("NKG7", "PRF1", "GZMB", "GNLY", "CTSW"),
  IFN_response = c("STAT1", "IRF1", "CXCL10", "ISG15", "IFIT3"),
  Antigen_presentation = c("HLA-DRA", "HLA-DRB1", "CD74", "B2M", "HLA-DPA1"),
  Exhaustion = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX")
)

hallmark <- list(
  HALLMARK_INTERFERON_GAMMA_RESPONSE = c("STAT1", "IRF1", "ISG15", "IFIT3", "CXCL10", "B2M"),
  HALLMARK_INFLAMMATORY_RESPONSE = c("IL1B", "TNF", "NFKB1", "CXCL8", "PTGS2", "CCL2"),
  HALLMARK_APOPTOSIS = c("BAX", "BCL2", "CASP3", "FAS", "TP53", "BID"),
  HALLMARK_OXIDATIVE_PHOSPHORYLATION = c("NDUFA1", "NDUFA2", "COX4I1", "ATP5F1A", "UQCRC1"),
  HALLMARK_IL6_JAK_STAT3_SIGNALING = c("IL6", "JAK1", "JAK2", "STAT3", "SOCS3")
)

# Shared gene universe
marker_genes <- unique(c(unlist(immune_small), unlist(hallmark)))
other_genes <- sprintf("GENE%04d", seq_len(900))
all_genes <- unique(c(marker_genes, other_genes))

# toy_expr: very small deterministic matrix + metadata
n_genes_toy <- 120
n_cells_toy <- 60
toy_genes <- all_genes[seq_len(n_genes_toy)]
toy_cells <- sprintf("toy_cell_%03d", seq_len(n_cells_toy))

toy_dense <- matrix(
  rpois(n_genes_toy * n_cells_toy, lambda = 2),
  nrow = n_genes_toy,
  ncol = n_cells_toy,
  dimnames = list(toy_genes, toy_cells)
)

# Inject group/celltype signal for a few pathways
toy_group <- rep(c("control", "treated"), each = n_cells_toy / 2)
toy_celltype <- sample(c("CD8_T", "B_cell", "Monocyte"), n_cells_toy, replace = TRUE)
toy_sample <- rep(sprintf("S%02d", 1:6), each = 10)

boost_genes <- intersect(c("PRF1", "GZMB", "NKG7"), rownames(toy_dense))
idx_boost <- which(toy_group == "treated" & toy_celltype == "CD8_T")
if (length(boost_genes) > 0 && length(idx_boost) > 0) {
  toy_dense[boost_genes, idx_boost] <- toy_dense[boost_genes, idx_boost] + 4
}

toy_expr <- list(
  expr = toy_dense,
  meta = data.frame(
    cell_id = toy_cells,
    sample = toy_sample,
    group = toy_group,
    celltype = toy_celltype,
    row.names = toy_cells,
    stringsAsFactors = FALSE
  )
)

# pbmc_small: sparse matrix + metadata
n_genes_pbmc <- 800
n_cells_pbmc <- 1000
pbmc_genes <- all_genes[seq_len(n_genes_pbmc)]
pbmc_cells <- sprintf("pbmc_cell_%04d", seq_len(n_cells_pbmc))

sp <- Matrix::rsparsematrix(nrow = n_genes_pbmc, ncol = n_cells_pbmc, density = 0.08)
sp@x <- pmax(0, round(sp@x * 3 + rnorm(length(sp@x), mean = 2, sd = 1.5), 0))
sp <- Matrix::drop0(sp)
rownames(sp) <- pbmc_genes
colnames(sp) <- pbmc_cells

sample_ids <- sprintf("PBMC_S%02d", rep(1:10, each = 100))
group_map <- setNames(rep(c("control", "treated"), each = 5), sprintf("PBMC_S%02d", 1:10))
pbmc_group <- unname(group_map[sample_ids])
pbmc_celltype <- sample(c("CD4_T", "CD8_T", "B_cell", "Monocyte", "NK"), n_cells_pbmc, replace = TRUE)

# Inject mild signal for treated CD8_T
sig_genes <- intersect(c("PRF1", "GZMB", "NKG7", "IFIT3", "STAT1"), rownames(sp))
idx_sig <- which(pbmc_group == "treated" & pbmc_celltype == "CD8_T")
if (length(sig_genes) > 0 && length(idx_sig) > 0) {
  m <- as.matrix(sp[sig_genes, idx_sig, drop = FALSE])
  m <- m + 1
  sp[sig_genes, idx_sig] <- Matrix::Matrix(m, sparse = TRUE)
}

pbmc_small <- list(
  expr = sp,
  meta = data.frame(
    cell_id = pbmc_cells,
    sample = sample_ids,
    group = pbmc_group,
    celltype = pbmc_celltype,
    row.names = pbmc_cells,
    stringsAsFactors = FALSE
  )
)

save(immune_small, file = "data/immune_small.rda", compress = "xz")
save(hallmark, file = "data/hallmark.rda", compress = "xz")
save(pbmc_small, file = "data/pbmc_small.rda", compress = "xz")
save(toy_expr, file = "data/toy_expr.rda", compress = "xz")

message("Data files written to data/*.rda")
