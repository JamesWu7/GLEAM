#' End-to-end scPathway workflow
#'
#' @param object Seurat object.
#' @param expr Expression matrix for matrix mode.
#' @param meta Metadata for matrix mode.
#' @param geneset Geneset input.
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
#' @param top_n Number of top pathways.
#' @param verbose Print messages.
#'
#' @return List containing score object and comparison result.
#' @export
run_scpathway <- function(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset = "hallmark",
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = "ensemble",
  group,
  sample = NULL,
  celltype = NULL,
  comparison = c("group", "celltype", "within_celltype"),
  target_celltype = NULL,
  level = "sample",
  top_n = 20,
  verbose = TRUE
) {
  comparison <- match.arg(comparison)

  sc <- score_pathway(
    object = object,
    expr = expr,
    meta = meta,
    geneset = geneset,
    seurat = seurat,
    assay = assay,
    layer = layer,
    slot = slot,
    method = method,
    verbose = verbose
  )

  tst <- switch(
    comparison,
    group = test_pathway(sc, group = group, sample = sample, level = level, verbose = verbose),
    celltype = compare_celltypes(sc, celltype = celltype, group = group, verbose = verbose),
    within_celltype = {
      if (is.null(target_celltype)) {
        stop("`target_celltype` is required when comparison = 'within_celltype'.", call. = FALSE)
      }
      compare_groups_within_celltype(
        score = sc,
        group = group,
        celltype = celltype,
        target_celltype = target_celltype,
        sample = sample,
        level = if (level %in% c("cell", "sample")) level else "sample",
        verbose = verbose
      )
    }
  )

  ord <- order(tst$table$p_adj, decreasing = FALSE, na.last = TRUE)
  top_tbl <- tst$table[head(ord, top_n), , drop = FALSE]

  list(score = sc, test = tst, top_table = top_tbl)
}
