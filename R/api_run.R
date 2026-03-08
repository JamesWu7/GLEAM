#' End-to-end GLEAM workflow
#'
#' @param object Seurat object.
#' @param expr Expression matrix for matrix mode.
#' @param meta Metadata for matrix mode.
#' @param geneset Geneset input.
#' @param geneset_source Geneset source.
#' @param seurat Input mode.
#' @param assay Seurat assay.
#' @param layer Seurat layer.
#' @param slot Legacy Seurat slot fallback.
#' @param method Scoring method.
#' @param group Group variable.
#' @param sample Sample variable.
#' @param celltype Celltype variable.
#' @param comparison Comparison mode.
#' @param target_celltype Target celltype for within-celltype comparison.
#' @param level Testing level.
#' @param pseudotime Pseudotime vector or metadata column name.
#' @param lineage Lineage vector or metadata column name.
#' @param region Spatial region/domain column name or vector.
#' @param backend Trajectory backend used when `comparison = 'trajectory'`.
#' @param top_n Number of top pathways.
#' @param verbose Print messages.
#'
#' @return List containing score object and comparison result.
#' @export
run_gleam <- function(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset = "hallmark",
  geneset_source = "auto",
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = "ensemble",
  group = NULL,
  sample = NULL,
  celltype = NULL,
  comparison = c("group", "celltype", "within_celltype", "trajectory", "spatial"),
  target_celltype = NULL,
  level = "sample",
  pseudotime = NULL,
  lineage = NULL,
  region = NULL,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot"),
  top_n = 20,
  verbose = TRUE
) {
  comparison <- match.arg(comparison)
  backend <- match.arg(backend)

  sc <- score_signature(
    object = object,
    expr = expr,
    meta = meta,
    geneset = geneset,
    geneset_source = geneset_source,
    seurat = seurat,
    assay = assay,
    layer = layer,
    slot = slot,
    method = method,
    verbose = verbose
  )

  tst <- switch(
    comparison,
    group = test_signature(sc, group = group, sample = sample, celltype = celltype, level = level, verbose = verbose),
    celltype = compare_celltypes(sc, celltype = celltype, group = group, verbose = verbose),
    within_celltype = {
      if (is.null(target_celltype)) stop("`target_celltype` is required when comparison = 'within_celltype'.", call. = FALSE)
      compare_groups_within_celltype(
        score = sc,
        group = group,
        celltype = celltype,
        target_celltype = target_celltype,
        sample = sample,
        level = if (level %in% c("cell", "sample")) level else "sample",
        verbose = verbose
      )
    },
    trajectory = {
      test_signature_trajectory(sc, pathway = NULL, pseudotime = pseudotime, lineage = lineage, backend = backend, verbose = verbose)
    },
    spatial = {
      if (is.null(region)) stop("`region` is required when comparison = 'spatial'.", call. = FALSE)
      test_signature(sc, group = group, sample = sample, region = region, level = "region", verbose = verbose)
    }
  )

  tbl <- if (inherits(tst, "gleam_test") || inherits(tst, "scpathway_test")) tst$table else tst
  ord <- order(tbl$p_adj, decreasing = FALSE, na.last = TRUE)
  top_tbl <- tbl[utils::head(ord, min(top_n, nrow(tbl))), , drop = FALSE]

  list(score = sc, test = tst, top_table = top_tbl)
}
