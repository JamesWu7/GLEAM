out_dir <- file.path("docs", "assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required to generate homepage figures.")
}

library(GLEAM)

data("toy_expr", package = "GLEAM")
data("spatial_medium_expr", package = "GLEAM")
data("spatial_medium_meta", package = "GLEAM")
data("spatial_medium_coords", package = "GLEAM")

sc <- score_signature(
  expr = toy_expr$expr,
  meta = toy_expr$meta,
  geneset = "immune_small",
  seurat = FALSE,
  method = "rank",
  min_genes = 3,
  verbose = FALSE
)

# Build a lightweight embedding from matrix data for figure generation.
pc <- stats::prcomp(t(as.matrix(toy_expr$expr)), center = TRUE, scale. = TRUE)
emb <- pc$x[, 1:2, drop = FALSE]
colnames(emb) <- c("UMAP_1", "UMAP_2")

p1 <- plot_embedding_score(sc, pathway = rownames(sc$score)[1], embedding = emb) +
  ggplot2::labs(title = "scRNA signature score on embedding")
ggplot2::ggsave(file.path(out_dir, "scrna_embedding_signature.png"), p1, width = 7.2, height = 5.2, dpi = 150)

p2 <- plot_dot_bar(sc, by = c("group", "celltype"), pathway = rownames(sc$score)[1:5]) +
  ggplot2::labs(title = "Dot-bar signature comparison")
ggplot2::ggsave(file.path(out_dir, "scrna_dotbar_signature.png"), p2, width = 8.0, height = 5.4, dpi = 150)

sp <- score_signature(
  expr = spatial_medium_expr,
  meta = spatial_medium_meta,
  geneset = "immune_small",
  seurat = FALSE,
  method = "rank",
  min_genes = 3,
  verbose = FALSE
)
p3 <- plot_spatial_score(sp, pathway = rownames(sp$score)[1], coords = spatial_medium_coords) +
  ggplot2::labs(title = "Spatial signature mapping")
ggplot2::ggsave(file.path(out_dir, "spatial_signature_map.png"), p3, width = 7.2, height = 5.2, dpi = 150)

p4 <- plot_pseudotime_score(sc, pathway = rownames(sc$score)[1], pseudotime = "pseudotime", lineage = "lineage") +
  ggplot2::labs(title = "Trajectory-aware signature trend")
ggplot2::ggsave(file.path(out_dir, "trajectory_signature_trend.png"), p4, width = 7.2, height = 5.2, dpi = 150)

message("[GLEAM] homepage figures generated in ", out_dir)
