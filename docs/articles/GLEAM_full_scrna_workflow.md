# GLEAM full scRNA-seq workflow

## 1) Build a Seurat object with standard preprocessing

For realistic scRNA analysis, signature scoring is most useful after a
standard cell-state workflow (normalization, variable feature selection,
dimensionality reduction, neighbor graph construction, and clustering).

``` r
library(Seurat)

ifnb_path <- system.file("extdata", "full_examples", "ifnb_seurat.rds", package = "GLEAM")
if (ifnb_path == "") ifnb_path <- file.path("inst", "extdata", "full_examples", "ifnb_seurat.rds")
seu <- readRDS(ifnb_path)
pick_first_col <- function(candidates, cols) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0L) return(NULL)
  hit[[1]]
}

stratified_keep <- function(meta, n_target, strata_candidates = character()) {
  if (nrow(meta) <= n_target) return(rownames(meta))
  strata <- intersect(strata_candidates, colnames(meta))
  all_ids <- rownames(meta)
  if (length(strata) == 0L) return(all_ids[seq_len(n_target)])

  key <- interaction(meta[, strata, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(all_ids, key)
  per_group <- max(1L, floor(n_target / max(1L, length(groups))))
  keep <- unlist(lapply(groups, function(ids) {
    ids <- sort(ids)
    head(ids, per_group)
  }), use.names = FALSE)
  if (length(keep) < n_target) {
    keep <- c(keep, head(setdiff(all_ids, keep), n_target - length(keep)))
  }
  unique(keep)[seq_len(min(n_target, length(unique(keep))))]
}

if (ncol(seu) > 7000) {
  md0 <- seu@meta.data
  keep <- stratified_keep(
    meta = md0,
    n_target = 7000,
    strata_candidates = c("orig.ident", "stim", "seurat_annotations", "seurat_clusters")
  )
  seu <- subset(seu, cells = keep)
}

md <- seu@meta.data
sample_col <- pick_first_col(c("sample", "orig.ident"), colnames(md))
if (is.null(sample_col)) {
  md$sample <- "sample_1"
  sample_col <- "sample"
} else if (sample_col != "sample") {
  md$sample <- as.character(md[[sample_col]])
}

group_col <- pick_first_col(c("stim", "group"), colnames(md))
if (is.null(group_col)) {
  md$group <- ifelse(seq_len(nrow(md)) <= nrow(md) / 2, "A", "B")
  group_col <- "group"
} else if (group_col != "group") {
  md$group <- as.character(md[[group_col]])
}

celltype_col <- pick_first_col(c("seurat_annotations", "celltype", "seurat_clusters"), colnames(md))
if (is.null(celltype_col)) {
  md$celltype <- as.character(Idents(seu))
  celltype_col <- "celltype"
} else if (celltype_col == "seurat_clusters") {
  md$celltype <- paste0("cluster_", md$seurat_clusters)
} else if (celltype_col != "celltype") {
  md$celltype <- as.character(md[[celltype_col]])
}
if (length(unique(md$sample)) < 2L) {
  md$sample <- ifelse(seq_len(nrow(md)) <= nrow(md) / 2, "sample_A", "sample_B")
}
if (length(unique(md$group)) < 2L) {
  s1 <- unique(md$sample)[1]
  md$group <- ifelse(md$sample == s1, "A", "B")
}
seu@meta.data <- md

# Basic QC context for tutorial reproducibility.
if (any(grepl("^MT-", rownames(seu)))) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
}
qc_cols <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(seu@meta.data))
summary(seu@meta.data[, qc_cols, drop = FALSE])

# Minimal, conventional Seurat preprocessing path.
if (!"pca" %in% names(seu@reductions)) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
}
ElbowPlot(seu, ndims = 30)
if (!"umap" %in% names(seu@reductions)) {
  seu <- FindNeighbors(seu, dims = 1:20)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = 1:20)
}
if (!"tsne" %in% names(seu@reductions)) {
  seu <- RunTSNE(seu, dims = 1:20)
}

if (!"pseudotime" %in% colnames(seu@meta.data)) {
  seu$pseudotime <- rank(Embeddings(seu, "pca")[, 1], ties.method = "average") / ncol(seu)
}
if (!"lineage" %in% colnames(seu@meta.data)) {
  seu$lineage <- as.character(seu[[celltype_col, drop = TRUE]])
}

# Show object context before downstream analysis
dim(seu)
seu
head(seu@meta.data[, unique(c("sample", group_col, celltype_col, "pseudotime", "lineage"))])
table(Idents(seu))
DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.35)
DimPlot(seu, reduction = "umap", group.by = "celltype", pt.size = 0.35)
DimPlot(seu, reduction = "tsne", group.by = "group", pt.size = 0.35)
```

## 2) Built-in geneset scoring + optional custom scoring

``` r
map_geneset_to_expr <- function(gs, expr_genes) {
  expr_genes <- as.character(expr_genes)
  expr_upper <- toupper(expr_genes)
  mapped <- lapply(gs, function(g) {
    idx <- match(toupper(unique(as.character(g))), expr_upper, nomatch = 0L)
    unique(expr_genes[idx[idx > 0L]])
  })
  mapped[vapply(mapped, length, integer(1)) > 0L]
}

fallback_signatures <- function(expr_genes) {
  g <- unique(as.character(expr_genes))
  g <- g[!is.na(g) & nzchar(g)]
  k <- min(30L, max(3L, floor(length(g) / 4)))
  idx1 <- seq_len(min(k, length(g)))
  idx2 <- seq.int(max(1L, length(g) - k + 1L), length(g))
  list(
    Signature_A = g[idx1],
    Signature_B = g[idx2]
  )
}

gs_h <- NULL
for (sp in c("human", "mouse")) {
  gs_try <- try(get_geneset("hallmark", source = "builtin", species = sp), silent = TRUE)
  if (inherits(gs_try, "try-error")) next
  gs_try <- map_geneset_to_expr(gs_try, rownames(seu))
  gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
  if (length(gs_try) > 0L) {
    gs_h <- gs_try
    break
  }
}
if (is.null(gs_h) || length(gs_h) == 0L) {
  message("[GLEAM] No matched Hallmark signatures; using object-derived fallback signatures.")
  gs_h <- fallback_signatures(rownames(seu))
}

sc <- score_signature(
  object = seu,
  geneset = gs_h,
  geneset_source = "list",
  seurat = TRUE,
  assay = "RNA",
  layer = "data",
  slot = "data",
  method = "ensemble",
  min_genes = 3
)

custom <- list(IFN_custom = c("STAT1", "IRF1", "ISG15", "IFIT3", "CXCL10"))
sc_custom <- score_signature(object = seu, geneset = custom, geneset_source = "list", seurat = TRUE, method = "mean", min_genes = 2)
```

## 3) Differential analysis after scoring

``` r
# cell-level exploratory
res_cell <- test_signature(sc, group = group_col, level = "cell", method = "wilcox")

# sample-level
res_sample <- test_signature(sc, group = group_col, sample = "sample", level = "sample", method = "wilcox")

# pseudobulk sample + celltype
res_pb <- test_signature(sc, group = group_col, sample = "sample", celltype = celltype_col, level = "pseudobulk", method = "wilcox")

# within-celltype (pick a celltype that has >=2 groups; otherwise skip gracefully)
ct_order <- names(sort(table(sc$meta[[celltype_col]]), decreasing = TRUE))
target_ct <- NA_character_
for (ct in ct_order) {
  idx <- sc$meta[[celltype_col]] == ct
  if (length(unique(sc$meta[[group_col]][idx])) >= 2L) {
    target_ct <- ct
    break
  }
}
if (!is.na(target_ct)) {
  wct_level <- if (length(unique(sc$meta$sample[sc$meta[[celltype_col]] == target_ct])) >= 2L) "sample" else "cell"
  res_wct <- test_signature(
    sc,
    group = group_col,
    sample = "sample",
    celltype = celltype_col,
    target_celltype = target_ct,
    level = wct_level
  )
} else {
  message("[GLEAM] within-celltype test skipped: no celltype has >=2 groups.")
  res_wct <- res_cell
}
```

## 4) Embedding and Seurat-style visualizations

``` r
top_pw <- res_pb$table$pathway[order(res_pb$table$p_adj)][1]
make_discrete_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(9, "Set1")
    if (n <= length(base)) return(base[seq_len(n)])
    return(grDevices::colorRampPalette(base)(n))
  }
  get_palette("gleam_discrete", n = n)
}

pal_group <- make_discrete_palette(length(unique(as.character(sc$meta[[group_col]]))))
pal_celltype <- make_discrete_palette(length(unique(as.character(sc$meta[[celltype_col]]))))

safe_plot <- function(label, expr) {
  p <- tryCatch(eval.parent(substitute(expr)), error = function(e) {
    message("[GLEAM] ", label, " skipped: ", conditionMessage(e))
    NULL
  })
  if (is.null(p)) return(invisible(NULL))
  tryCatch(print(p), error = function(e) {
    message("[GLEAM] ", label, " render skipped: ", conditionMessage(e))
  })
  invisible(NULL)
}

sig_col <- paste0("GLEAM_signature_", gsub("[^A-Za-z0-9_]+", "_", top_pw))
sig_vals <- rep(NA_real_, ncol(seu))
names(sig_vals) <- colnames(seu)
sig_cells <- intersect(colnames(seu), colnames(sc$score))
sig_vals[sig_cells] <- as.numeric(sc$score[top_pw, sig_cells])
seu[[sig_col]] <- sig_vals

safe_plot(
  "seurat_featureplot",
  Seurat::FeaturePlot(
    seu,
    features = sig_col,
    reduction = "umap",
    cols = get_palette("gleam_continuous", n = 128, continuous = TRUE),
    order = TRUE,
    pt.size = 0.8
  ) + ggplot2::labs(title = "Seurat FeaturePlot-style signature map")
)
safe_plot(
  "seurat_vlnplot",
  Seurat::VlnPlot(
    seu,
    features = sig_col,
    group.by = celltype_col,
    split.by = group_col,
    pt.size = 0
  ) + Seurat::NoLegend() + ggplot2::labs(title = "Seurat VlnPlot-style signature distribution")
)

safe_plot("embedding_umap", plot_embedding_score(sc, pathway = top_pw, object = seu, reduction = "umap"))
safe_plot("embedding_pca", plot_embedding_score(sc, pathway = top_pw, object = seu, reduction = "pca"))
safe_plot("embedding_tsne", plot_embedding_score(sc, pathway = top_pw, object = seu, reduction = "tsne"))

safe_plot("dot", plot_dot(sc, by = c("group", "celltype")))
safe_plot("dot_bar", plot_dot_bar(sc, by = c("group", "celltype"), pathway = top_pw, color_palette = "gleam_continuous"))
safe_plot("violin", plot_violin(sc, pathway = rownames(sc$score)[1], group = group_col, palette = pal_group))
safe_plot("split_violin", plot_split_violin(sc, pathway = rownames(sc$score)[1], x = celltype_col, split.by = group_col, palette = pal_group))
safe_plot("ridge", plot_ridge(sc, pathway = rownames(sc$score)[1], group = celltype_col, palette = pal_celltype))
safe_plot("pseudobulk_box", plot_pseudobulk_box(sc, pathway = rownames(sc$score)[1], group = group_col, sample = "sample", celltype = celltype_col, palette = pal_group))
safe_plot("volcano", plot_volcano(res_pb))
```

## 5) Trajectory / pseudotime integration

``` r
# Generic metadata-driven trajectory fields
traj_res <- test_signature_trajectory(sc, pseudotime = "pseudotime", lineage = "lineage", method = "spearman")

if (!exists("make_discrete_palette")) {
  make_discrete_palette <- function(n) {
    n <- max(1L, as.integer(n))
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      base <- RColorBrewer::brewer.pal(9, "Set1")
      if (n <= length(base)) return(base[seq_len(n)])
      return(grDevices::colorRampPalette(base)(n))
    }
    get_palette("gleam_discrete", n = n)
  }
}

if (!exists("safe_plot")) {
  safe_plot <- function(label, expr) {
    p <- tryCatch(eval.parent(substitute(expr)), error = function(e) {
      message("[GLEAM] ", label, " skipped: ", conditionMessage(e))
      NULL
    })
    if (is.null(p)) return(invisible(NULL))
    tryCatch(print(p), error = function(e) {
      message("[GLEAM] ", label, " render skipped: ", conditionMessage(e))
    })
    invisible(NULL)
  }
}

pal_lineage <- make_discrete_palette(length(unique(as.character(sc$meta$lineage))))
umap_emb <- Seurat::Embeddings(seu, "umap")

safe_plot(
  "pseudotime",
  plot_pseudotime_score(
    sc,
    pathway = rownames(sc$score)[1],
    pseudotime = "pseudotime",
    lineage = "lineage",
    palette = pal_lineage
  )
)
safe_plot(
  "trajectory_umap",
  plot_trajectory_score(
    sc,
    pathway = rownames(sc$score)[1],
    embeddings = umap_emb,
    reduction = "umap"
  )
)
```

### Optional Monocle2 / Monocle3 / slingshot object pathways

``` r
# Monocle2 object (package monocle)
# traj_res_m2 <- test_signature_trajectory(sc, pseudotime = cds_monocle2, lineage = cds_monocle2, method = "spearman")

# Monocle3 object (package monocle3)
# traj_res_m3 <- test_signature_trajectory(sc, pseudotime = cds_monocle3, lineage = "lineage", method = "spearman")

# slingshot object
# traj_res_sl <- test_signature_trajectory(sc, pseudotime = sds, lineage = sds, method = "spearman")
```

## 6) Export score tables and compare methods

``` r
long_tbl <- pivot_scores_long(sc)
sum_tbl <- summarize_scores(sc, by = c("sample", "group"), fun = "mean")
out_csv <- file.path(tempdir(), "gleam_scores_scrna.csv")
export_scores(sc, out_csv, format = "csv", include_meta = TRUE)

head(long_tbl)
head(sum_tbl)
```
