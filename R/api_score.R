#' Score pathway activity per cell
#'
#' Supports Seurat v5 input (`seurat = TRUE`) and direct matrix/sparse matrix
#' input (`seurat = FALSE`). Seurat is optional and only required for Seurat mode.
#'
#' @param object Seurat object when `seurat = TRUE`.
#' @param expr Expression matrix or `dgCMatrix` when `seurat = FALSE`.
#' @param meta Metadata data.frame for matrix mode.
#' @param geneset Geneset input: built-in name, named list, GMT path, or data.frame.
#' @param seurat Logical input mode.
#' @param assay Seurat assay name.
#' @param layer Seurat v5 layer name. If `NULL`, tries `data` then `counts`.
#' @param slot Legacy Seurat slot fallback.
#' @param method Scoring method.
#' @param min_genes Minimum genes per pathway.
#' @param max_genes Maximum genes per pathway.
#' @param auc_max_rank AUC top-rank proportion.
#' @param ensemble_methods Methods used by ensemble.
#' @param ensemble_combine Ensemble combine strategy.
#' @param verbose Whether to print messages.
#'
#' @return An object of class `scpathway_score`.
#' @export
score_pathway <- function(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset,
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = c("rank", "auc", "zmean", "ensemble"),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "auc", "zmean"),
  ensemble_combine = c("mean", "median"),
  verbose = TRUE
) {
  method <- match.arg(method)
  ensemble_combine <- match.arg(ensemble_combine)

  check_input_mode(object = object, expr = expr, seurat = seurat)

  expr_use <- extract_expr(
    object = object,
    expr = expr,
    assay = assay,
    layer = layer,
    slot = slot,
    seurat = seurat
  )

  meta_use <- extract_meta(
    object = object,
    meta = meta,
    expr = expr_use,
    seurat = seurat
  )

  check_expr_meta(expr_use, meta_use)
  if (!identical(rownames(meta_use), colnames(expr_use))) {
    meta_use <- meta_use[colnames(expr_use), , drop = FALSE]
  }

  gs <- get_geneset(geneset)
  gs <- check_geneset(gs, min_genes = min_genes, max_genes = max_genes)
  matched <- match_geneset(gs, expr_genes = rownames(expr_use), verbose = verbose)
  gs_match <- matched$geneset

  keep <- vapply(gs_match, length, integer(1)) > 0L
  if (!any(keep)) {
    stop("No geneset overlaps with expression matrix genes.", call. = FALSE)
  }
  gs_match <- gs_match[keep]

  score_mat <- switch(
    method,
    rank = score_rank_matrix(expr_use, gs_match, verbose = verbose),
    auc = score_auc_matrix(expr_use, gs_match, auc_max_rank = auc_max_rank, verbose = verbose),
    zmean = score_zmean_matrix(expr_use, gs_match, verbose = verbose),
    ensemble = score_ensemble_matrix(
      expr_use,
      gs_match,
      methods = ensemble_methods,
      combine = ensemble_combine,
      auc_max_rank = auc_max_rank,
      verbose = verbose
    )
  )

  new_scpathway_score(
    score = score_mat,
    meta = meta_use,
    method = method,
    geneset_name = if (is.character(geneset) && length(geneset) == 1L) geneset else "custom",
    geneset_info = matched$info,
    params = list(
      seurat = seurat,
      assay = assay,
      layer = layer,
      slot = slot,
      min_genes = min_genes,
      max_genes = max_genes,
      auc_max_rank = auc_max_rank,
      ensemble_methods = ensemble_methods,
      ensemble_combine = ensemble_combine
    )
  )
}
