out_dir <- file.path("assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (requireNamespace("GLEAM", quietly = TRUE)) {
  library(GLEAM)
} else if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", export_all = FALSE, quiet = TRUE)
} else {
  stop("Either installed package 'GLEAM' or package 'pkgload' is required.")
}

generate_fallback <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to generate homepage figures.")
  }

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
  pc <- stats::prcomp(t(as.matrix(toy_expr$expr)), center = TRUE, scale. = TRUE)
  emb <- pc$x[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1", "UMAP_2")

  sp <- score_signature(
    expr = spatial_medium_expr,
    meta = spatial_medium_meta,
    geneset = "immune_small",
    seurat = FALSE,
    method = "rank",
    min_genes = 3,
    verbose = FALSE
  )

  tissue_bg <- as.raster(matrix(colorRampPalette(c("#f8f5ea", "#e8dcc4", "#d2b48c"))(256), nrow = 16))

  p1 <- plot_embedding_score(sc, pathway = rownames(sc$score)[1], embedding = emb) +
    ggplot2::labs(title = "Signature score on embedding")
  ggplot2::ggsave(file.path(out_dir, "embedding_signature_feature.png"), p1, width = 8.0, height = 5.6, dpi = 160)

  p2 <- plot_spatial_score(sp, pathway = rownames(sp$score)[1], coords = spatial_medium_coords, image = tissue_bg) +
    ggplot2::labs(title = "Signature score on spatial slice")
  ggplot2::ggsave(file.path(out_dir, "spatial_slice_signature.png"), p2, width = 8.0, height = 5.6, dpi = 160)

  p3 <- plot_dot_bar(sc, by = c("group", "celltype"), pathway = rownames(sc$score)[1:5]) +
    ggplot2::labs(title = "Dot-bar signature comparison")
  ggplot2::ggsave(file.path(out_dir, "signature_dotbar_compare.png"), p3, width = 8.0, height = 5.4, dpi = 160)

  p4 <- plot_pseudotime_score(sc, pathway = rownames(sc$score)[1], pseudotime = "pseudotime", lineage = "lineage") +
    ggplot2::labs(title = "Trajectory-aware signature trend")
  ggplot2::ggsave(file.path(out_dir, "trajectory_signature_trend.png"), p4, width = 7.8, height = 5.2, dpi = 160)
}

can_render_vignette <- requireNamespace("rmarkdown", quietly = TRUE) &&
  requireNamespace("knitr", quietly = TRUE) &&
  isTRUE(rmarkdown::pandoc_available())

if (can_render_vignette) {
  vignette_file <- file.path("vignettes", "GLEAM_homepage_showcase.Rmd")
  if (!file.exists(vignette_file)) stop("Missing vignette file: ", vignette_file)

  tmp_dir <- file.path(tempdir(), "gleam_homepage_figures")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  ok <- TRUE
  tryCatch(
    rmarkdown::render(
      input = vignette_file,
      output_dir = tmp_dir,
      quiet = TRUE,
      knit_root_dir = normalizePath("."),
      envir = new.env(parent = globalenv())
    ),
    error = function(e) {
      ok <<- FALSE
      message("[GLEAM] vignette render failed, falling back to direct workflow generation: ", conditionMessage(e))
    }
  )

  if (ok) {
    fig_dir <- file.path(tmp_dir, "GLEAM_homepage_showcase_files", "figure-html")
    if (dir.exists(fig_dir)) {
      copy_named <- function(src_name, dest_name) {
        src <- file.path(fig_dir, src_name)
        if (!file.exists(src)) stop("Expected vignette figure not found: ", src)
        file.copy(src, file.path(out_dir, dest_name), overwrite = TRUE)
      }
      copy_named("embedding_feature-1.png", "embedding_signature_feature.png")
      copy_named("spatial_slice-1.png", "spatial_slice_signature.png")
      copy_named("dotbar_compare-1.png", "signature_dotbar_compare.png")
      copy_named("trajectory_trend-1.png", "trajectory_signature_trend.png")
      message("[GLEAM] homepage figures generated from vignette output in ", out_dir)
    } else {
      message("[GLEAM] vignette figure directory missing, using fallback generation.")
      generate_fallback()
    }
  } else {
    generate_fallback()
  }
} else {
  message("[GLEAM] rmarkdown/knitr/pandoc unavailable; using fallback generation.")
  generate_fallback()
}
