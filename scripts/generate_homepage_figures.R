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
  "trajectory_signature_trend.png"
)

has_existing_figures <- function() {
  all(file.exists(file.path(out_dir, figure_targets)))
}

map_geneset_to_expr <- function(gs, expr_genes) {
  expr_genes <- as.character(expr_genes)
  expr_upper <- toupper(expr_genes)
  mapped <- lapply(gs, function(g) {
    idx <- match(toupper(unique(as.character(g))), expr_upper, nomatch = 0L)
    unique(expr_genes[idx[idx > 0L]])
  })
  mapped[vapply(mapped, length, integer(1)) > 0L]
}

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

make_discrete_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(9, "Set1")
    if (n <= length(base)) return(base[seq_len(n)])
    return(grDevices::colorRampPalette(base)(n))
  }
  get_palette("gleam_discrete", n = n)
}

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
    keep <- stratified_keep(
      meta = md0,
      n_target = 5000,
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
    gs_try <- map_geneset_to_expr(gs_try, rownames(seu))
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      hallmark_gs <- gs_try
      break
    }
  }
  if (is.null(hallmark_gs) || length(hallmark_gs) == 0L) {
    hallmark_gs <- fallback_signatures(rownames(seu))
  }

  sc <- score_signature(
    object = seu,
    geneset = hallmark_gs,
    geneset_source = "list",
    seurat = TRUE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  sp_obj <- readRDS(a1_path)
  if (ncol(sp_obj) > 3500) {
    keep <- stratified_keep(
      meta = sp_obj@meta.data,
      n_target = 3500,
      strata_candidates = c("orig.ident", "sample", "seurat_annotations", "seurat_clusters")
    )
    sp_obj <- subset(sp_obj, cells = keep)
  }
  md_sp <- sp_obj@meta.data
  region_col <- pick_first_col(c("seurat_annotations", "region", "seurat_clusters"), colnames(md_sp))
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
    gs_try <- map_geneset_to_expr(gs_try, rownames(sp_obj))
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      gs_sp <- gs_try
      break
    }
  }
  if (is.null(gs_sp) || length(gs_sp) == 0L) {
    gs_sp <- fallback_signatures(rownames(sp_obj))
  }
  sp <- score_signature(
    object = sp_obj,
    geneset = gs_sp,
    geneset_source = "list",
    seurat = TRUE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  p1 <- plot_embedding_score(sc, pathway = rownames(sc$score)[1], object = seu, reduction = "umap") +
    ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 10.0, height = 7.0, dpi = 220)

  p2 <- plot_spatial_score(sp, pathway = rownames(sp$score)[1], object = sp_obj) +
    ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 11.0, height = 8.2, dpi = 240)

  p3 <- plot_dot_bar(sc, by = c(group_col, celltype_col), pathway = rownames(sc$score)[1:5]) +
    ggplot2::labs(title = "Dot-bar signature comparison")
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 12.8, height = 8.8, dpi = 220)

  pal_lineage <- make_discrete_palette(length(unique(as.character(sc$meta$lineage))))
  p4 <- tryCatch(
    plot_pseudotime_score(
      sc,
      pathway = rownames(sc$score)[1],
      pseudotime = "pseudotime",
      lineage = "lineage",
      palette = pal_lineage
    ) + ggplot2::labs(title = "Trajectory-aware signature trend"),
    error = function(e) {
      message("[GLEAM] trajectory trend fallback: ", conditionMessage(e))
      plot_pseudotime_score(
        sc,
        pathway = rownames(sc$score)[1],
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
    gs_try <- map_geneset_to_expr(gs_try, rownames(expr_pbmc))
    gs_try <- gs_try[vapply(gs_try, length, integer(1)) >= 3L]
    if (length(gs_try) > 0L) {
      gs_pbmc <- gs_try
      break
    }
  }
  if (is.null(gs_pbmc) || length(gs_pbmc) == 0L) {
    gs_pbmc <- fallback_signatures(rownames(expr_pbmc))
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
  ct_top <- head(ct_top, 5L)
  keep_cells <- rownames(meta_pbmc)[meta_pbmc$celltype %in% ct_top]
  sc_vis <- sc
  sc_vis$score <- sc$score[, keep_cells, drop = FALSE]
  sc_vis$meta <- sc$meta[keep_cells, , drop = FALSE]

  gs_sp <- fallback_signatures(rownames(expr_sp))
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

  p1 <- plot_embedding_score(sc, pathway = rownames(sc$score)[1], embedding = emb, reduction = "umap") +
    ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 10.0, height = 7.0, dpi = 220)

  p2 <- plot_spatial_score(sc_sp, pathway = rownames(sc_sp$score)[1], coords = coords_sp, image = tissue_bg) +
    ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 11.0, height = 8.2, dpi = 240)

  sig_n <- min(3L, nrow(sc_vis$score))
  sc_dot <- sc_vis
  sc_dot$score <- sc_vis$score[seq_len(sig_n), , drop = FALSE]
  rownames(sc_dot$score) <- paste0("Sig_", seq_len(sig_n))
  p3 <- plot_dot_bar(sc_dot, by = c("group", "celltype"))
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 12.0, height = 7.8, dpi = 220)

  n_lin <- length(unique(as.character(sc_vis$meta$lineage)))
  pal_lineage <- make_discrete_palette(n_lin)
  p4 <- plot_pseudotime_score(
    sc_vis,
    pathway = rownames(sc_vis$score)[1],
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
  ok <- isTRUE(generate_from_builtin_examples())
}

if (!ok && has_existing_figures()) {
  message("[GLEAM] keeping existing homepage figures in ", out_dir)
  ok <- TRUE
}

if (!ok) {
  stop("Failed to generate homepage figures from full examples, built-in examples, and existing files.")
}
