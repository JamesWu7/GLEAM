out_dir <- file.path("man", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
} else if (requireNamespace("GLEAM", quietly = TRUE)) {
  library(GLEAM)
} else {
  stop("Either installed package 'GLEAM' or package 'pkgload' is required.")
}

figure_targets <- c(
  "embedding_signature_feature.png",
  "spatial_slice_signature.png",
  "signature_dotbar_compare.png",
  "signature_violin.png",
  "signature_split_violin.png",
  "signature_ridge.png",
  "trajectory_signature_trend.png"
)

has_existing_figures <- function() {
  all(file.exists(file.path(out_dir, figure_targets)))
}

.pick_first_col <- GLEAM:::.pick_first_col
.stratified_keep <- GLEAM:::.stratified_keep
.map_geneset_to_expr <- GLEAM:::.map_geneset_to_expr
.fallback_signatures <- GLEAM:::.fallback_signatures
.make_discrete_palette <- GLEAM:::.make_discrete_palette

find_full_example <- function(file_name) {
  p <- system.file("extdata", "full_examples", file_name, package = "GLEAM")
  if (p == "") p <- file.path("inst", "extdata", "full_examples", file_name)
  if (!file.exists(p)) stop("Missing full example file: ", file_name)
  p
}

copy_vignette_figures <- function(fig_dir) {
  copy_named <- function(src_name, dest_name) {
    src <- file.path(fig_dir, src_name)
    if (!file.exists(src)) stop("Expected vignette figure not found: ", src)
    file.copy(src, file.path(out_dir, dest_name), overwrite = TRUE)
  }

  copy_named("embedding_feature-1.png", "embedding_signature_feature.png")
  copy_named("spatial_slice-1.png", "spatial_slice_signature.png")
  copy_named("dotbar_compare-1.png", "signature_dotbar_compare.png")
  copy_named("violin_compare-1.png", "signature_violin.png")
  copy_named("split_violin_compare-1.png", "signature_split_violin.png")
  if (file.exists(file.path(fig_dir, "ridge_compare-1.png"))) {
    copy_named("ridge_compare-1.png", "signature_ridge.png")
  } else {
    copy_named("violin_compare-1.png", "signature_ridge.png")
  }
  copy_named("trajectory_trend-1.png", "trajectory_signature_trend.png")
}

generate_from_full_examples <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    message("[GLEAM] Seurat/ggplot2 unavailable; cannot regenerate homepage figures from full examples.")
    return(FALSE)
  }

  library(Seurat)

  ifnb_path <- find_full_example("ifnb_seurat.rds")
  a1_path <- find_full_example("stxBrain_anterior1_seurat.rds")

  seu <- readRDS(ifnb_path)
  md0 <- seu@meta.data
  if (ncol(seu) > 5000) {
    keep <- .stratified_keep(
      meta = md0,
      n_target = 5000,
      strata_candidates = c("orig.ident", "stim", "seurat_annotations", "seurat_clusters")
    )
    seu <- subset(seu, cells = keep)
  }
  md <- seu@meta.data
  sample_col <- .pick_first_col(c("sample", "orig.ident"), colnames(md))
  if (is.null(sample_col)) {
    md$sample <- "sample_1"
    sample_col <- "sample"
  } else if (sample_col != "sample") {
    md$sample <- as.character(md[[sample_col]])
  }

  group_col <- .pick_first_col(c("stim", "group"), colnames(md))
  if (is.null(group_col)) {
    md$group <- ifelse(seq_len(nrow(md)) <= nrow(md) / 2, "A", "B")
    group_col <- "group"
  } else if (group_col != "group") {
    md$group <- as.character(md[[group_col]])
  }

  celltype_col <- .pick_first_col(c("seurat_annotations", "celltype", "seurat_clusters"), colnames(md))
  if (is.null(celltype_col)) {
    md$celltype <- as.character(Idents(seu))
    celltype_col <- "celltype"
  } else if (celltype_col == "seurat_clusters") {
    md$celltype <- paste0("cluster_", md$seurat_clusters)
  } else if (celltype_col != "celltype") {
    md$celltype <- as.character(md[[celltype_col]])
  }
  seu@meta.data <- md

  if (!"pca" %in% names(seu@reductions)) {
    seu <- NormalizeData(seu)
    seu <- FindVariableFeatures(seu)
    seu <- ScaleData(seu)
    seu <- RunPCA(seu)
  }
  if (!"umap" %in% names(seu@reductions)) {
    seu <- FindNeighbors(seu, dims = 1:20)
    seu <- FindClusters(seu, resolution = 0.4)
    seu <- RunUMAP(seu, dims = 1:20)
  }
  md <- seu@meta.data
  md$pseudotime <- rank(Embeddings(seu, "pca")[, 1], ties.method = "average") / ncol(seu)
  md$lineage <- as.character(md[[celltype_col]])
  seu@meta.data <- md

  hallmark_gs <- NULL
  for (sp in c("human", "mouse")) {
    gs_try <- try(get_geneset("hallmark", source = "builtin", species = sp), silent = TRUE)
    if (inherits(gs_try, "try-error")) next
    gs_try <- .map_geneset_to_expr(gs_try, rownames(seu), min_genes = 1L)
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      hallmark_gs <- gs_try
      break
    }
  }
  if (is.null(hallmark_gs) || length(hallmark_gs) == 0L) {
    hallmark_gs <- .fallback_signatures(rownames(seu), max_genes = 30L, min_genes = 3L)
  }

  sc <- score_signature(
    object = seu,
    geneset = hallmark_gs,
    geneset_source = "list",
    seurat = TRUE,
    layer = "counts",
    slot = "counts",
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  sp_obj <- readRDS(a1_path)
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    sp_obj <- tryCatch(SeuratObject::UpdateSeuratObject(sp_obj), error = function(e) sp_obj)
  }
  if (ncol(sp_obj) > 3500) {
    keep <- .stratified_keep(
      meta = sp_obj@meta.data,
      n_target = 3500,
      strata_candidates = c("orig.ident", "sample", "seurat_annotations", "seurat_clusters")
    )
    sp_obj <- subset(sp_obj, cells = keep)
  }
  md_sp <- sp_obj@meta.data
  region_col <- .pick_first_col(c("seurat_annotations", "region", "seurat_clusters"), colnames(md_sp))
  if (is.null(region_col)) {
    md_sp$region <- as.character(Idents(sp_obj))
  } else if (region_col == "seurat_clusters") {
    md_sp$region <- paste0("cluster_", md_sp$seurat_clusters)
  } else if (region_col != "region") {
    md_sp$region <- as.character(md_sp[[region_col]])
  }
  sp_obj@meta.data <- md_sp

  gs_sp <- NULL
  for (sp_name in c("mouse", "human")) {
    gs_try <- try(get_geneset("hallmark", source = "builtin", species = sp_name), silent = TRUE)
    if (inherits(gs_try, "try-error")) next
    gs_try <- .map_geneset_to_expr(gs_try, rownames(sp_obj), min_genes = 1L)
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      gs_sp <- gs_try
      break
    }
  }
  if (is.null(gs_sp) || length(gs_sp) == 0L) {
    gs_sp <- .fallback_signatures(rownames(sp_obj), max_genes = 30L, min_genes = 3L)
  }
  sp <- score_signature(
    object = sp_obj,
    geneset = gs_sp,
    geneset_source = "list",
    seurat = TRUE,
    layer = "counts",
    slot = "counts",
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  top_sig <- rownames(sc$score)[1]
  p1 <- plot_embedding_score(
    score = sc,
    pathway = top_sig,
    object = seu,
    reduction = "umap",
    point_size = 0.62
  ) + ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 10.0, height = 7.0, dpi = 220)

  p2 <- plot_spatial_score(
    score = sp,
    pathway = rownames(sp$score)[1],
    object = sp_obj,
    point_size = 1.65
  ) + ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 11.0, height = 8.2, dpi = 240)

  ct_rank <- names(sort(table(sc$meta[[celltype_col]]), decreasing = TRUE))
  ct_keep <- head(ct_rank, 8L)
  keep_cells <- rownames(sc$meta)[sc$meta[[celltype_col]] %in% ct_keep]
  sc_vis <- sc
  sc_vis$score <- sc$score[, keep_cells, drop = FALSE]
  sc_vis$meta <- sc$meta[keep_cells, , drop = FALSE]
  grp_levels <- unique(as.character(sc_vis$meta[[group_col]]))
  sc_vis$meta[[celltype_col]] <- factor(as.character(sc_vis$meta[[celltype_col]]), levels = ct_keep)
  sc_vis$meta[[group_col]] <- factor(as.character(sc_vis$meta[[group_col]]), levels = grp_levels)
  pal_celltype_vis <- setNames(get_palette("gleam_discrete", n = length(ct_keep), continuous = FALSE), ct_keep)
  pal_group_vis <- setNames(get_palette("gleam_discrete", n = length(grp_levels), continuous = FALSE), grp_levels)

  p3 <- plot_dot_bar(sc_vis, by = c(group_col, celltype_col), pathway = top_sig) +
    ggplot2::labs(title = "Dot-bar signature comparison")
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 12.8, height = 7.8, dpi = 220)

  p5 <- plot_violin(
    score = sc_vis,
    pathway = top_sig,
    group = celltype_col,
    palette = pal_celltype_vis,
    point_size = 0.25,
    alpha = 0.7
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(title = "Violin plot across major cell types")
  ggplot2::ggsave(file.path(out_dir, "signature_violin.png"), p5, width = 11.0, height = 6.4, dpi = 220)

  p6 <- plot_split_violin(
    score = sc_vis,
    pathway = top_sig,
    x = celltype_col,
    split.by = group_col,
    palette = pal_group_vis,
    alpha = 0.74
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(title = "Split violin by group within major cell types")
  ggplot2::ggsave(file.path(out_dir, "signature_split_violin.png"), p6, width = 11.4, height = 6.6, dpi = 220)

  p7 <- tryCatch(
    plot_ridge(sc_vis, pathway = top_sig, group = celltype_col, alpha = 0.72, palette = pal_celltype_vis) +
      ggplot2::labs(title = "Ridge distribution across major cell types"),
    error = function(e) {
      message("[GLEAM] ridge fallback to violin: ", conditionMessage(e))
      p5
    }
  )
  ggplot2::ggsave(file.path(out_dir, "signature_ridge.png"), p7, width = 11.0, height = 6.6, dpi = 220)

  pal_lineage <- .make_discrete_palette(length(unique(as.character(sc$meta$lineage))))
  p4 <- tryCatch(
    plot_pseudotime_score(
      sc_vis,
      pathway = top_sig,
      pseudotime = "pseudotime",
      lineage = "lineage",
      palette = pal_lineage
    ) + ggplot2::labs(title = "Trajectory-aware signature trend"),
    error = function(e) {
      message("[GLEAM] trajectory trend fallback: ", conditionMessage(e))
      plot_pseudotime_score(
        sc_vis,
        pathway = top_sig,
        pseudotime = "pseudotime",
        lineage = NULL,
        palette = pal_lineage
      ) + ggplot2::labs(title = "Trajectory-aware signature trend")
    }
  )
  ggplot2::ggsave(file.path(out_dir, "trajectory_signature_trend.png"), p4, width = 10.0, height = 6.8, dpi = 220)

  TRUE
}

load_example_data <- function(name) {
  env <- new.env(parent = emptyenv())
  data_file <- file.path("data", paste0(name, ".rda"))
  if (file.exists(data_file)) {
    load(data_file, envir = env)
  } else {
    data(list = name, package = "GLEAM", envir = env)
  }
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop("Missing example dataset: ", name)
  }
  get(name, envir = env, inherits = FALSE)
}

generate_from_builtin_examples <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("[GLEAM] ggplot2 unavailable; cannot generate homepage figures from built-in examples.")
    return(FALSE)
  }

  expr_pbmc <- load_example_data("pbmc_medium_matrix")
  meta_pbmc <- load_example_data("pbmc_medium_meta")
  expr_sp <- load_example_data("spatial_medium_expr")
  meta_sp <- load_example_data("spatial_medium_meta")
  coords_sp <- load_example_data("spatial_medium_coords")

  meta_pbmc <- as.data.frame(meta_pbmc, stringsAsFactors = FALSE)
  if ("cell_id" %in% colnames(meta_pbmc)) {
    rownames(meta_pbmc) <- meta_pbmc$cell_id
  }
  meta_pbmc <- meta_pbmc[colnames(expr_pbmc), , drop = FALSE]

  meta_sp <- as.data.frame(meta_sp, stringsAsFactors = FALSE)
  if ("cell_id" %in% colnames(meta_sp)) {
    rownames(meta_sp) <- meta_sp$cell_id
  }
  coords_sp <- as.data.frame(coords_sp, stringsAsFactors = FALSE)
  rownames(coords_sp) <- rownames(meta_sp)
  coords_sp <- coords_sp[colnames(expr_sp), c("x", "y"), drop = FALSE]
  meta_sp <- meta_sp[colnames(expr_sp), , drop = FALSE]

  gs_pbmc <- NULL
  for (sp_name in c("human", "mouse")) {
    gs_try <- try(get_geneset("hallmark", source = "builtin", species = sp_name), silent = TRUE)
    if (inherits(gs_try, "try-error")) next
    gs_try <- .map_geneset_to_expr(gs_try, rownames(expr_pbmc), min_genes = 1L)
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      gs_pbmc <- gs_try
      break
    }
  }
  if (is.null(gs_pbmc) || length(gs_pbmc) == 0L) {
    gs_pbmc <- .fallback_signatures(rownames(expr_pbmc), max_genes = 30L, min_genes = 3L)
  }

  sc <- score_signature(
    expr = expr_pbmc,
    meta = meta_pbmc,
    geneset = gs_pbmc,
    geneset_source = "list",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  top_idx <- head(order(Matrix::rowSums(expr_pbmc > 0), decreasing = TRUE), 220L)
  emb <- stats::prcomp(
    t(as.matrix(expr_pbmc[top_idx, , drop = FALSE])),
    rank. = 2,
    center = TRUE,
    scale. = TRUE
  )$x[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  rownames(emb) <- colnames(expr_pbmc)

  ct_top <- names(sort(table(meta_pbmc$celltype), decreasing = TRUE))
  ct_top <- head(ct_top, 8L)
  keep_cells <- rownames(meta_pbmc)[meta_pbmc$celltype %in% ct_top]
  sc_vis <- sc
  sc_vis$score <- sc$score[, keep_cells, drop = FALSE]
  sc_vis$meta <- sc$meta[keep_cells, , drop = FALSE]
  grp_levels <- unique(as.character(sc_vis$meta$group))
  sc_vis$meta$celltype <- factor(as.character(sc_vis$meta$celltype), levels = ct_top)
  sc_vis$meta$group <- factor(as.character(sc_vis$meta$group), levels = grp_levels)
  pal_celltype_vis <- setNames(get_palette("gleam_discrete", n = length(ct_top), continuous = FALSE), ct_top)
  pal_group_vis <- setNames(get_palette("gleam_discrete", n = length(grp_levels), continuous = FALSE), grp_levels)
  top_sig <- rownames(sc_vis$score)[1]

  gs_sp <- .fallback_signatures(rownames(expr_sp), max_genes = 30L, min_genes = 3L)
  sc_sp <- score_signature(
    expr = expr_sp,
    meta = meta_sp,
    geneset = gs_sp,
    geneset_source = "list",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  tissue_bg <- as.raster(matrix(
    colorRampPalette(c("#f7f4ea", "#eadfc8", "#d9c6a0", "#c8a97b"))(256),
    nrow = 16
  ))

  p1 <- plot_embedding_score(sc, pathway = top_sig, embedding = emb, reduction = "umap", point_size = 0.62) +
    ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 10.0, height = 7.0, dpi = 220)

  p2 <- plot_spatial_score(sc_sp, pathway = rownames(sc_sp$score)[1], coords = coords_sp, image = tissue_bg, point_size = 1.65) +
    ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 11.0, height = 8.2, dpi = 240)

  sig_n <- min(1L, nrow(sc_vis$score))
  sc_dot <- sc_vis
  sc_dot$score <- sc_vis$score[seq_len(sig_n), , drop = FALSE]
  rownames(sc_dot$score) <- paste0("Sig_", seq_len(sig_n))
  p3 <- plot_dot_bar(sc_dot, by = c("group", "celltype"))
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 12.0, height = 7.8, dpi = 220)

  p5 <- plot_violin(
    score = sc_vis,
    pathway = top_sig,
    group = "celltype",
    palette = pal_celltype_vis,
    point_size = 0.25,
    alpha = 0.7
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(title = "Violin plot across major cell types")
  ggplot2::ggsave(file.path(out_dir, "signature_violin.png"), p5, width = 11.0, height = 6.4, dpi = 220)

  p6 <- plot_split_violin(
    score = sc_vis,
    pathway = top_sig,
    x = "celltype",
    split.by = "group",
    palette = pal_group_vis,
    alpha = 0.74
  ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
    ggplot2::labs(title = "Split violin by group within major cell types")
  ggplot2::ggsave(file.path(out_dir, "signature_split_violin.png"), p6, width = 11.4, height = 6.6, dpi = 220)

  p7 <- tryCatch(
    plot_ridge(sc_vis, pathway = top_sig, group = "celltype", alpha = 0.72, palette = pal_celltype_vis) +
      ggplot2::labs(title = "Ridge distribution across major cell types"),
    error = function(e) {
      message("[GLEAM] ridge fallback to violin: ", conditionMessage(e))
      p5
    }
  )
  ggplot2::ggsave(file.path(out_dir, "signature_ridge.png"), p7, width = 11.0, height = 6.6, dpi = 220)

  n_lin <- length(unique(as.character(sc_vis$meta$lineage)))
  pal_lineage <- .make_discrete_palette(n_lin)
  p4 <- plot_pseudotime_score(
    sc_vis,
    pathway = top_sig,
    pseudotime = "pseudotime",
    lineage = "lineage",
    palette = pal_lineage
  ) + ggplot2::labs(title = "Trajectory-aware signature trend")
  ggplot2::ggsave(file.path(out_dir, "trajectory_signature_trend.png"), p4, width = 10.0, height = 6.8, dpi = 220)

  TRUE
}

can_render_vignette <- requireNamespace("rmarkdown", quietly = TRUE) &&
  requireNamespace("knitr", quietly = TRUE) &&
  isTRUE(rmarkdown::pandoc_available())

ok <- FALSE

if (can_render_vignette) {
  vignette_file <- file.path("scripts", "GLEAM_homepage_showcase.Rmd")
  if (!file.exists(vignette_file)) stop("Missing vignette file: ", vignette_file)

  tmp_dir <- file.path(tempdir(), "gleam_homepage_figures")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch(
    rmarkdown::render(
      input = vignette_file,
      output_dir = tmp_dir,
      quiet = TRUE,
      knit_root_dir = normalizePath("."),
      envir = new.env(parent = globalenv())
    ),
    error = function(e) {
      message("[GLEAM] vignette render failed, falling back to direct full-example generation: ", conditionMessage(e))
    }
  )

  fig_dir <- file.path(tmp_dir, "GLEAM_homepage_showcase_files", "figure-html")
  if (dir.exists(fig_dir)) {
    tryCatch(
      {
        copy_vignette_figures(fig_dir)
        message("[GLEAM] homepage figures generated from vignette output in ", out_dir)
        ok <- TRUE
      },
      error = function(e) {
        message("[GLEAM] failed to copy vignette figures, using direct full-example generation: ", conditionMessage(e))
      }
    )
  }
}

if (!ok) {
  ok <- isTRUE(generate_from_full_examples())
}

if (!ok) {
  if (!requireNamespace("Seurat", quietly = TRUE) && has_existing_figures()) {
    message("[GLEAM] Seurat unavailable; keeping existing homepage figures in ", out_dir)
    ok <- TRUE
  }
}

if (!ok) {
  ok <- isTRUE(generate_from_builtin_examples())
}

if (!ok && has_existing_figures()) {
  message("[GLEAM] keeping existing homepage figures in ", out_dir)
  ok <- TRUE
}

if (!ok) {
  stop("Failed to generate homepage figures from full examples, built-in examples, and existing files.")
}
