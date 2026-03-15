#' Score pathway activity per cell
#'
#' Supports Seurat v4/v5 input (`seurat = TRUE`) and direct matrix/sparse matrix
#' input (`seurat = FALSE`). Seurat is optional and only required for Seurat mode.
#'
#' @param object Seurat object when `seurat = TRUE`.
#' @param expr Expression matrix or `dgCMatrix` when `seurat = FALSE`.
#' @param meta Metadata data.frame for matrix mode.
#' @param geneset Geneset input. Also supports signed signatures:
#'   `list(signatureA = list(up = c(...), down = c(...)))`.
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
#' @param method Scoring method. Canonical methods include:
#'   `rank`, `mean`, `zscore`, `scaled_mean`, `robust_mean`, `ensemble`,
#'   `AddModuleScore`, `UCell`, `AUCell`, `ssGSEA`, `GSVA`, `singscore`.
#' @param min_genes Minimum genes per pathway.
#' @param max_genes Maximum genes per pathway.
#' @param auc_max_rank AUC top-rank proportion used by `AUCell`.
#' @param ensemble_methods Methods used by ensemble.
#' @param ensemble_combine Ensemble combine strategy.
#' @param ensemble_standardize Harmonization before ensemble aggregation:
#'   `zscore` (recommended) or `rank`.
#' @param ensemble_weights Optional named numeric vector of per-method weights.
#' @param method_params Optional named list for method-specific parameters.
#'   Supported keys include:
#'   `auc_max_rank` (AUCell), `alpha` or `ssgsea_alpha` (ssGSEA),
#'   `nbin`/`ctrl`/`seed` (AddModuleScore), `kcdf` (GSVA),
#'   `ucell_max_rank` (UCell).
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
    "rank", "mean", "zscore", "scaled_mean", "robust_mean", "ensemble",
    "AddModuleScore", "UCell", "AUCell", "ssGSEA", "GSVA", "singscore"
  ),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "zscore", "mean"),
  ensemble_combine = c("mean", "median"),
  ensemble_standardize = c("zscore", "rank"),
  ensemble_weights = NULL,
  method_params = list(),
  verbose = TRUE
) {
  geneset_source <- match.arg(geneset_source)
  ontology <- match.arg(ontology)
  ensemble_combine <- match.arg(ensemble_combine)
  ensemble_standardize <- match.arg(ensemble_standardize)
  method <- canonicalize_scoring_method(as.character(method)[1])
  ensemble_methods <- canonicalize_scoring_methods(ensemble_methods)

  if (!is.list(method_params)) {
    stop("`method_params` must be a named list.", call. = FALSE)
  }

  check_input_mode(object = object, expr = expr, seurat = seurat)

  expr_use <- extract_expr(
    object = object,
    expr = expr,
    assay = assay,
    layer = layer,
    slot = slot,
    seurat = seurat
  )
  expr_info <- attr(expr_use, "gleam_expr_info")

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
  gs_norm <- .normalize_signed_genesets(gs)
  base_meta <- get_geneset_metadata(gs)

  gs_combined <- lapply(gs_norm, function(x) unique(c(x$up, x$down)))
  gs_combined <- check_geneset(gs_combined, min_genes = min_genes, max_genes = max_genes)
  pathways <- names(gs_combined)
  gs_norm <- gs_norm[pathways]
  gs_up <- lapply(gs_norm, `[[`, "up")
  gs_down <- lapply(gs_norm, `[[`, "down")

  expr_genes <- rownames(expr_use)
  gs_up_match <- lapply(gs_up, intersect, y = expr_genes)
  gs_down_match <- lapply(gs_down, intersect, y = expr_genes)
  n_input <- vapply(gs_combined, length, integer(1))
  n_up_input <- vapply(gs_up, length, integer(1))
  n_down_input <- vapply(gs_down, length, integer(1))
  n_up_match <- vapply(gs_up_match, length, integer(1))
  n_down_match <- vapply(gs_down_match, length, integer(1))
  n_match <- n_up_match + n_down_match

  keep <- n_match > 0L
  if (!any(keep)) {
    stop("No geneset overlaps with expression matrix genes.", call. = FALSE)
  }

  pathways <- pathways[keep]
  gs_up_match <- gs_up_match[pathways]
  gs_down_match <- gs_down_match[pathways]
  n_input <- n_input[pathways]
  n_match <- n_match[pathways]
  n_up_input <- n_up_input[pathways]
  n_down_input <- n_down_input[pathways]
  n_up_match <- n_up_match[pathways]
  n_down_match <- n_down_match[pathways]

  genes_dropped <- lapply(pathways, function(pw) {
    orig <- unique(c(gs_up[[pw]], gs_down[[pw]]))
    hit <- unique(c(gs_up_match[[pw]], gs_down_match[[pw]]))
    setdiff(orig, hit)
  })
  names(genes_dropped) <- pathways

  info <- data.frame(
    pathway = pathways,
    n_genes_input = as.integer(n_input),
    n_genes_matched = as.integer(n_match),
    frac_matched = as.numeric(n_match / pmax(1L, n_input)),
    n_genes_up_input = as.integer(n_up_input),
    n_genes_down_input = as.integer(n_down_input),
    n_genes_up_matched = as.integer(n_up_match),
    n_genes_down_matched = as.integer(n_down_match),
    stringsAsFactors = FALSE
  )
  if (!is.null(base_meta) && "pathway" %in% colnames(base_meta)) {
    info <- merge(base_meta, info, by = "pathway", all.y = TRUE, sort = FALSE)
    info <- info[match(pathways, info$pathway), , drop = FALSE]
  }

  if (verbose) {
    message(sprintf("[GLEAM] matched pathways: %d", length(pathways)))
    message(sprintf("[GLEAM] median matched genes: %.1f", stats::median(n_match)))
  }

  .score_unsigned <- function(meth, genesets_use, params_use) {
    if (length(genesets_use) == 0L) {
      return(matrix(numeric(0), nrow = 0, ncol = ncol(expr_use)))
    }
    if (meth != "ensemble") check_scoring_method(meth)

    if (meth == "rank") {
      return(score_rank_matrix(expr_use, genesets_use, verbose = verbose))
    }
    if (meth == "zscore") {
      return(score_zscore_matrix(expr_use, genesets_use, verbose = verbose))
    }
    if (meth == "mean") {
      return(score_mean_matrix(expr_use, genesets_use, verbose = verbose))
    }
    if (meth == "scaled_mean") {
      return(score_scaled_mean_matrix(expr_use, genesets_use, verbose = verbose))
    }
    if (meth == "robust_mean") {
      return(score_robust_matrix(expr_use, genesets_use, verbose = verbose))
    }

    score_optional_wrapper(
      method = meth,
      object = object,
      expr = expr_use,
      genesets = genesets_use,
      seurat = seurat,
      assay = assay,
      slot = slot %||% "data",
      method_params = params_use
    )
  }

  .score_signed <- function(meth) {
    params_use <- .resolve_method_params(method = meth, method_params = method_params, auc_max_rank = auc_max_rank)
    out <- matrix(0, nrow = length(pathways), ncol = ncol(expr_use), dimnames = list(pathways, colnames(expr_use)))

    up_nonempty <- pathways[vapply(gs_up_match[pathways], length, integer(1)) > 0L]
    down_nonempty <- pathways[vapply(gs_down_match[pathways], length, integer(1)) > 0L]

    if (length(up_nonempty) > 0L) {
      up_mat <- .score_unsigned(meth, gs_up_match[up_nonempty], params_use)
      out[rownames(up_mat), ] <- out[rownames(up_mat), , drop = FALSE] + up_mat
    }
    if (length(down_nonempty) > 0L) {
      down_mat <- .score_unsigned(meth, gs_down_match[down_nonempty], params_use)
      out[rownames(down_mat), ] <- out[rownames(down_mat), , drop = FALSE] - down_mat
    }

    list(score = out, params = params_use)
  }

  method_parameters_used <- NULL
  score_mat <- if (method == "ensemble") {
    em <- setdiff(ensemble_methods, "ensemble")
    if (length(em) == 0L) {
      stop("`ensemble_methods` must include at least one non-ensemble method.", call. = FALSE)
    }
    scored <- lapply(em, .score_signed)
    mats <- lapply(scored, `[[`, "score")
    names(mats) <- em
    params_by_method <- lapply(scored, `[[`, "params")
    names(params_by_method) <- em
    method_parameters_used <- list(
      ensemble_methods = em,
      ensemble_combine = ensemble_combine,
      ensemble_standardize = ensemble_standardize,
      ensemble_weights = ensemble_weights,
      method_params_by_method = params_by_method
    )

    mats_h <- lapply(mats, function(x) .harmonize_score_matrix(x, mode = ensemble_standardize))
    w <- .resolve_ensemble_weights(methods = names(mats_h), weights = ensemble_weights)
    arr <- simplify2array(mats_h)

    out <- if (ensemble_combine == "mean" || !is.null(ensemble_weights)) {
      apply(arr, c(1, 2), function(v) {
        ok <- !is.na(v)
        if (!any(ok)) return(NA_real_)
        ww <- w[ok]
        if (all(ww == 0)) return(NA_real_)
        sum(v[ok] * ww) / sum(ww)
      })
    } else {
      apply(arr, c(1, 2), function(v) stats::median(v, na.rm = TRUE))
    }
    rownames(out) <- rownames(mats_h[[1]])
    colnames(out) <- colnames(mats_h[[1]])
    out
  } else {
    scored <- .score_signed(method)
    method_parameters_used <- scored$params
    scored$score
  }

  new_gleam_score(
    score = score_mat,
    meta = meta_use,
    method = method,
    geneset_name = if (is.character(geneset) && length(geneset) == 1L) geneset else "custom",
    geneset_info = info,
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
      ensemble_standardize = ensemble_standardize,
      ensemble_weights = ensemble_weights,
      method_params = method_params,
      geneset_size_original = n_input,
      geneset_size_matched = n_match,
      genes_dropped = genes_dropped,
      expression_layer_used = expr_info,
      method_parameters_used = method_parameters_used
    )
  )
}

#' @keywords internal
.normalize_signed_genesets <- function(gs) {
  out <- lapply(gs, function(x) {
    up <- character(0)
    down <- character(0)

    if (is.list(x) && !is.data.frame(x)) {
      nm <- tolower(names(x) %||% rep("", length(x)))
      if (any(nm == "up") || any(nm == "down")) {
        if (any(nm == "up")) up <- x[[which(nm == "up")[1]]]
        if (any(nm == "down")) down <- x[[which(nm == "down")[1]]]
      } else {
        up <- unlist(x, use.names = FALSE)
      }
    } else {
      up <- x
    }

    up <- unique(as.character(up))
    down <- unique(as.character(down))
    up <- up[!is.na(up) & nzchar(up)]
    down <- down[!is.na(down) & nzchar(down)]
    list(up = up, down = down)
  })

  names(out) <- names(gs)
  out
}

#' @keywords internal
.resolve_method_params <- function(method, method_params, auc_max_rank) {
  mp <- method_params %||% list()
  method_key <- tolower(method)
  nested <- mp[[method]]
  if (is.null(nested)) nested <- mp[[method_key]]
  if (is.list(nested)) {
    mp <- utils::modifyList(mp, nested)
  }

  if (method == "AUCell") {
    return(list(auc_max_rank = mp$auc_max_rank %||% auc_max_rank))
  }
  if (method == "ssGSEA") {
    return(list(
      alpha = mp$alpha %||% mp$ssgsea_alpha %||% 0.25,
      normalize = mp$normalize %||% TRUE
    ))
  }
  if (method == "AddModuleScore") {
    return(list(
      nbin = mp$nbin %||% 24,
      ctrl = mp$ctrl %||% 100,
      seed = mp$seed %||% 1
    ))
  }
  if (method == "GSVA") {
    return(list(kcdf = mp$kcdf %||% "Gaussian"))
  }
  if (method == "UCell") {
    return(list(ucell_max_rank = mp$ucell_max_rank %||% NULL))
  }
  list()
}
