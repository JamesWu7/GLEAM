set.seed(20260307)
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required.")

# Base genesets ------------------------------------------------------------
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

marker_genes <- unique(c(unlist(immune_small), unlist(hallmark)))
other_genes <- sprintf("GENE%05d", seq_len(4500))
all_genes <- unique(c(marker_genes, other_genes))

# toy_expr -----------------------------------------------------------------
n_genes_toy <- 120
n_cells_toy <- 60
toy_genes <- all_genes[seq_len(n_genes_toy)]
toy_cells <- sprintf("toy_cell_%03d", seq_len(n_cells_toy))

toy_dense <- matrix(rpois(n_genes_toy * n_cells_toy, lambda = 2), nrow = n_genes_toy, ncol = n_cells_toy,
                    dimnames = list(toy_genes, toy_cells))
toy_group <- rep(c("control", "treated"), each = n_cells_toy / 2)
toy_celltype <- sample(c("CD8_T", "B_cell", "Monocyte"), n_cells_toy, replace = TRUE)
toy_sample <- rep(sprintf("S%02d", 1:6), each = 10)
toy_pseudotime <- seq(0, 1, length.out = n_cells_toy)

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
    donor = rep(c("D1", "D2", "D3"), each = 20),
    batch = rep(c("B1", "B2"), each = 30),
    condition = toy_group,
    pseudotime = toy_pseudotime,
    lineage = ifelse(toy_celltype %in% c("CD8_T", "B_cell"), "L1", "L2"),
    row.names = toy_cells,
    stringsAsFactors = FALSE
  )
)

# pbmc_small (kept for backward examples) ----------------------------------
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

# pbmc_medium_matrix/meta ---------------------------------------------------
n_genes_med <- 2200
n_cells_med <- 3200
med_genes <- all_genes[seq_len(n_genes_med)]
med_cells <- sprintf("pbmc_med_cell_%05d", seq_len(n_cells_med))

pbmc_medium_matrix <- Matrix::rsparsematrix(nrow = n_genes_med, ncol = n_cells_med, density = 0.06)
pbmc_medium_matrix@x <- pmax(0, round(pbmc_medium_matrix@x * 4 + rnorm(length(pbmc_medium_matrix@x), mean = 2.5, sd = 1.2), 0))
pbmc_medium_matrix <- Matrix::drop0(pbmc_medium_matrix)
rownames(pbmc_medium_matrix) <- med_genes
colnames(pbmc_medium_matrix) <- med_cells

sample_med <- sprintf("M_S%02d", rep(1:16, each = n_cells_med / 16))
group_med_map <- setNames(rep(c("control", "treated"), each = 8), sprintf("M_S%02d", 1:16))
group_med <- unname(group_med_map[sample_med])
celltype_med <- sample(c("CD4_T", "CD8_T", "B_cell", "Monocyte", "NK", "DC"), n_cells_med, replace = TRUE,
                       prob = c(0.25, 0.18, 0.18, 0.2, 0.15, 0.04))
donor_med <- sprintf("Donor_%02d", rep(1:8, each = n_cells_med / 8))
batch_med <- rep(c("Batch1", "Batch2", "Batch3", "Batch4"), each = n_cells_med / 4)
condition_med <- ifelse(group_med == "treated", "stim", "unstim")
pseudotime_med <- pmin(1, pmax(0, stats::runif(n_cells_med, min = 0, max = 1) + ifelse(celltype_med %in% c("CD8_T", "NK"), 0.15, -0.05)))
lineage_med <- ifelse(celltype_med %in% c("CD4_T", "CD8_T", "NK"), "Lymphoid", "Myeloid")

# add signal
signal_genes <- intersect(c("PRF1", "GZMB", "NKG7", "STAT1", "IFIT3", "ISG15"), rownames(pbmc_medium_matrix))
idx_sig <- which(group_med == "treated" & celltype_med %in% c("CD8_T", "NK"))
if (length(signal_genes) > 0 && length(idx_sig) > 0) {
  m <- as.matrix(pbmc_medium_matrix[signal_genes, idx_sig, drop = FALSE])
  m <- m + 2
  pbmc_medium_matrix[signal_genes, idx_sig] <- Matrix::Matrix(m, sparse = TRUE)
}

pbmc_medium_meta <- data.frame(
  cell_id = med_cells,
  sample = sample_med,
  group = group_med,
  celltype = celltype_med,
  donor = donor_med,
  batch = batch_med,
  condition = condition_med,
  pseudotime = pseudotime_med,
  lineage = lineage_med,
  row.names = med_cells,
  stringsAsFactors = FALSE
)

# spatial example -----------------------------------------------------------
n_spots <- 1800
n_genes_sp <- 1500
sp_genes <- all_genes[seq_len(n_genes_sp)]
sp_cells <- sprintf("spatial_spot_%04d", seq_len(n_spots))

spatial_medium_expr <- Matrix::rsparsematrix(nrow = n_genes_sp, ncol = n_spots, density = 0.07)
spatial_medium_expr@x <- pmax(0, round(spatial_medium_expr@x * 3 + rnorm(length(spatial_medium_expr@x), mean = 1.8, sd = 1.1), 0))
spatial_medium_expr <- Matrix::drop0(spatial_medium_expr)
rownames(spatial_medium_expr) <- sp_genes
colnames(spatial_medium_expr) <- sp_cells

grid_n <- ceiling(sqrt(n_spots))
coords <- expand.grid(x = seq_len(grid_n), y = seq_len(grid_n))
coords <- coords[seq_len(n_spots), , drop = FALSE]
rownames(coords) <- sp_cells

section <- rep(c("Section1", "Section2"), each = n_spots / 2)
region <- ifelse(coords$x < median(coords$x) & coords$y < median(coords$y), "R1", ifelse(coords$x > median(coords$x), "R2", "R3"))
sample_sp <- sprintf("SP_S%02d", rep(1:12, each = n_spots / 12))
group_sp <- ifelse(sample_sp %in% sprintf("SP_S%02d", 1:6), "control", "treated")
domain <- ifelse(region == "R1", "immune_edge", ifelse(region == "R2", "stromal_core", "mixed_zone"))

spatial_medium_meta <- data.frame(
  cell_id = sp_cells,
  sample = sample_sp,
  group = group_sp,
  section = section,
  region = region,
  domain = domain,
  x = coords$x,
  y = coords$y,
  pseudotime = pmin(1, pmax(0, coords$x / max(coords$x) + stats::rnorm(n_spots, 0, 0.05))),
  lineage = ifelse(region == "R1", "Lymphoid", "Myeloid"),
  row.names = sp_cells,
  stringsAsFactors = FALSE
)
spatial_medium_coords <- coords

# Save ----------------------------------------------------------------------
if (!dir.exists("data")) dir.create("data", recursive = TRUE)
save(immune_small, file = "data/immune_small.rda", compress = "xz")
save(hallmark, file = "data/hallmark.rda", compress = "xz")
save(toy_expr, file = "data/toy_expr.rda", compress = "xz")
save(pbmc_small, file = "data/pbmc_small.rda", compress = "xz")
save(pbmc_medium_matrix, file = "data/pbmc_medium_matrix.rda", compress = "xz")
save(pbmc_medium_meta, file = "data/pbmc_medium_meta.rda", compress = "xz")
save(spatial_medium_expr, file = "data/spatial_medium_expr.rda", compress = "xz")
save(spatial_medium_meta, file = "data/spatial_medium_meta.rda", compress = "xz")
save(spatial_medium_coords, file = "data/spatial_medium_coords.rda", compress = "xz")

message("[GLEAM] data/*.rda refreshed")
