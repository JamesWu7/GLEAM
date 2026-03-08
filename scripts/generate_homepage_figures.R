out_dir <- file.path("assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("GLEAM", quietly = TRUE)) {
  library(GLEAM)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
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
  p1_path <- find_full_example("stxBrain_posterior1_seurat.rds")

  seu <- readRDS(ifnb_path)
  if (ncol(seu) > 5000) seu <- subset(seu, cells = colnames(seu)[1:5000])
  md <- seu@meta.data
  if (!"sample" %in% colnames(md)) md$sample <- if ("orig.ident" %in% colnames(md)) as.character(md$orig.ident) else "sample_1"
  if (!"group" %in% colnames(md)) md$group <- if ("stim" %in% colnames(md)) as.character(md$stim) else ifelse(seq_len(nrow(md)) <= nrow(md) / 2, "A", "B")
  if (!"celltype" %in% colnames(md)) md$celltype <- if ("seurat_clusters" %in% colnames(md)) paste0("cluster_", md$seurat_clusters) else as.character(Idents(seu))
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
  md$lineage <- md$celltype
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

  sp_a1 <- readRDS(a1_path)
  sp_p1 <- readRDS(p1_path)
  sp_seu <- merge(sp_a1, y = sp_p1, add.cell.ids = c("anterior", "posterior"))
  if (ncol(sp_seu) > 6000) sp_seu <- subset(sp_seu, cells = colnames(sp_seu)[1:6000])
  sp_md <- sp_seu@meta.data
  if (!"sample" %in% colnames(sp_md)) sp_md$sample <- if ("orig.ident" %in% colnames(sp_md)) as.character(sp_md$orig.ident) else "sample_1"
  if (!"group" %in% colnames(sp_md)) sp_md$group <- ifelse(grepl("anterior", sp_md$sample, ignore.case = TRUE), "anterior", "posterior")
  if (!"region" %in% colnames(sp_md)) sp_md$region <- if ("seurat_clusters" %in% colnames(sp_md)) paste0("cluster_", sp_md$seurat_clusters) else sp_md$group
  if (!"pca" %in% names(sp_seu@reductions)) {
    sp_seu <- NormalizeData(sp_seu)
    sp_seu <- FindVariableFeatures(sp_seu)
    sp_seu <- ScaleData(sp_seu)
    sp_seu <- RunPCA(sp_seu)
  }
  sp_seu@meta.data <- sp_md

  sp_genes <- rownames(sp_seu)
  half_n <- max(30L, floor(length(sp_genes) / 2))
  idx_a <- seq_len(min(half_n, length(sp_genes)))
  idx_b <- seq(from = min(half_n + 1L, length(sp_genes)), to = length(sp_genes))
  gs_sp <- list(
    Spatial_signature_A = unique(sp_genes[idx_a])[1:min(30L, length(unique(sp_genes[idx_a])))],
    Spatial_signature_B = unique(sp_genes[idx_b])[1:min(30L, length(unique(sp_genes[idx_b])))]
  )

  sp <- score_signature(
    object = sp_seu,
    geneset = gs_sp,
    geneset_source = "list",
    seurat = TRUE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  coords <- if (all(c("x", "y") %in% colnames(sp_md))) {
    data.frame(x = sp_md$x, y = sp_md$y, row.names = rownames(sp_md))
  } else if (all(c("imagecol", "imagerow") %in% colnames(sp_md))) {
    data.frame(x = sp_md$imagecol, y = sp_md$imagerow, row.names = rownames(sp_md))
  } else {
    emb <- Embeddings(sp_seu, reduction = "pca")[, 1:2, drop = FALSE]
    data.frame(x = emb[, 1], y = emb[, 2], row.names = rownames(emb))
  }

  tissue_bg <- as.raster(matrix(colorRampPalette(c("#f8f5ea", "#e8dcc4", "#d2b48c"))(256), nrow = 16))

  p1 <- plot_embedding_score(sc, pathway = rownames(sc$score)[1], object = seu, reduction = "umap") +
    ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 8.0, height = 5.6, dpi = 160)

  p2 <- plot_spatial_score(sp, pathway = rownames(sp$score)[1], coords = coords, image = tissue_bg) +
    ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 8.0, height = 5.6, dpi = 160)

  p3 <- plot_dot_bar(sc, by = c("group", "celltype"), pathway = rownames(sc$score)[1:5]) +
    ggplot2::labs(title = "Dot-bar signature comparison")
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 8.0, height = 5.4, dpi = 160)

  p4 <- plot_pseudotime_score(sc, pathway = rownames(sc$score)[1], pseudotime = "pseudotime", lineage = "lineage") +
    ggplot2::labs(title = "Trajectory-aware signature trend")
  ggplot2::ggsave(file.path(out_dir, "trajectory_signature_trend.png"), p4, width = 7.8, height = 5.2, dpi = 160)

  TRUE
}

can_render_vignette <- requireNamespace("rmarkdown", quietly = TRUE) &&
  requireNamespace("knitr", quietly = TRUE) &&
  isTRUE(rmarkdown::pandoc_available())

ok <- FALSE

if (can_render_vignette) {
  vignette_file <- file.path("vignettes", "GLEAM_homepage_showcase.Rmd")
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

if (!ok && has_existing_figures()) {
  message("[GLEAM] keeping existing homepage figures in ", out_dir)
  ok <- TRUE
}

if (!ok) {
  stop("Failed to generate homepage figures. Install Seurat or provide pre-generated files under assets/figures.")
}
