#' Score pathway activity per cell
#'
#' Supports Seurat v4/v5 input (`seurat = TRUE`) and direct matrix/sparse matrix
#' input (`seurat = FALSE`). Seurat is optional and only required for Seurat mode.
#'
#' @param object Seurat object when `seurat = TRUE`.
#' @param expr Expression matrix or `dgCMatrix` when `seurat = FALSE`.
#' @param meta Metadata data.frame for matrix mode.
#' @param geneset Geneset input.
#' @param geneset_source Geneset source (`auto`, `builtin`, `list`, `gmt`, `data.frame`,
#'   `msigdb`, `go`, `kegg`, `reactome`).
#' @param species Species label for source-aware geneset loading.
#' @param collection Collection parameter for source-aware geneset loading.
#' @param subcollection Subcollection parameter for source-aware geneset loading.
#' @param ontology GO ontology (`BP`, `MF`, `CC`) for GO source.
#' @param seurat Logical input mode.
#' @param assay Seurat assay name.
#' @param layer Seurat v5 layer name. If `NULL`, tries `data` then `counts`.
#' @param slot Legacy Seurat slot fallback.
#' @param method Scoring method. See [list_scoring_methods()].
#' @param min_genes Minimum genes per pathway.
#' @param max_genes Maximum genes per pathway.
#' @param auc_max_rank AUC top-rank proportion.
#' @param ensemble_methods Methods used by ensemble.
#' @param ensemble_combine Ensemble combine strategy.
#' @param method_params Optional named list for method-specific parameters.
#' @param verbose Whether to print messages.
#'
#' @return An object of class `gleam_score`.
#' @keywords internal
#' @noRd
score_pathway <- function(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset,
  geneset_source = c("auto", "builtin", "list", "gmt", "data.frame", "msigdb", "go", "kegg", "reactome"),
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  ontology = c("BP", "MF", "CC"),
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = c(
    "rank", "auc", "zmean", "mean", "scaled_mean", "robust", "singscore_like", "ssgsea_like", "ensemble",
    "addmodulescore", "ucell", "aucell", "gsva", "singscore"
  ),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "auc", "zmean"),
  ensemble_combine = c("mean", "median"),
  method_params = list(),
  verbose = TRUE
) {
  geneset_source <- match.arg(geneset_source)
  ontology <- match.arg(ontology)
  method <- tolower(match.arg(method))
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

  gs <- get_geneset(
    geneset = geneset,
    source = geneset_source,
    species = species,
    collection = collection,
    subcollection = subcollection,
    ontology = ontology
  )
  gs <- check_geneset(gs, min_genes = min_genes, max_genes = max_genes)
  matched <- match_geneset(gs, expr_genes = rownames(expr_use), verbose = verbose)
  gs_match <- matched$geneset

  keep <- vapply(gs_match, length, integer(1)) > 0L
  if (!any(keep)) {
    stop("No geneset overlaps with expression matrix genes.", call. = FALSE)
  }
  gs_match <- gs_match[keep]

  .score_one <- function(meth) {
    meth <- tolower(meth)
    if (meth != "ensemble") check_scoring_method(meth)

    if (meth == "rank") {
      return(score_rank_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "auc") {
      max_rank <- method_params$auc_max_rank %||% auc_max_rank
      return(score_auc_matrix(expr_use, gs_match, auc_max_rank = max_rank, verbose = verbose))
    }
    if (meth == "zmean") {
      return(score_zmean_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "mean") {
      return(score_mean_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "scaled_mean") {
      return(score_scaled_mean_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "robust") {
      return(score_robust_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "singscore_like") {
      return(score_singscore_like_matrix(expr_use, gs_match, verbose = verbose))
    }
    if (meth == "ssgsea_like") {
      alpha <- method_params$ssgsea_alpha %||% 0.25
      return(score_ssgsea_like_matrix(expr_use, gs_match, alpha = alpha, verbose = verbose))
    }

    score_optional_wrapper(
      method = meth,
      object = object,
      expr = expr_use,
      genesets = gs_match,
      seurat = seurat,
      assay = assay,
      slot = slot %||% "data"
    )
  }

  score_mat <- if (method == "ensemble") {
    mats <- lapply(ensemble_methods, .score_one)
    arr <- simplify2array(mats)
    out <- if (ensemble_combine == "mean") {
      apply(arr, c(1, 2), function(v) mean(v, na.rm = TRUE))
    } else {
      apply(arr, c(1, 2), function(v) stats::median(v, na.rm = TRUE))
    }
    rownames(out) <- rownames(mats[[1]])
    colnames(out) <- colnames(mats[[1]])
    out
  } else {
    .score_one(method)
  }

  new_scpathway_score(
    score = score_mat,
    meta = meta_use,
    method = method,
    geneset_name = if (is.character(geneset) && length(geneset) == 1L) geneset else "custom",
    geneset_info = matched$info,
    params = list(
      seurat = seurat,
      seurat_version = if (seurat) detect_seurat_version(object) else NA_integer_,
      assay = assay,
      layer = layer,
      slot = slot,
      geneset_source = geneset_source,
      species = species,
      collection = collection,
      subcollection = subcollection,
      ontology = ontology,
      min_genes = min_genes,
      max_genes = max_genes,
      auc_max_rank = auc_max_rank,
      ensemble_methods = ensemble_methods,
      ensemble_combine = ensemble_combine,
      method_params = method_params
    )
  )
}
