#' @keywords internal
score_optional_wrapper <- function(method, object = NULL, expr = NULL, genesets, seurat = FALSE, assay = NULL, slot = "data", method_params = list()) {
  method <- canonicalize_scoring_method(method)

  if (method == "AddModuleScore") {
    check_method_dependency(method, "Seurat")
    if (!seurat || is.null(object)) {
      stop("`AddModuleScore` requires `seurat = TRUE` and a Seurat object.", call. = FALSE)
    }
    nbin <- method_params$nbin %||% 24
    ctrl <- method_params$ctrl %||% 100
    seed <- method_params$seed %||% 1
    out <- Seurat::AddModuleScore(
      object = object,
      features = unname(genesets),
      assay = assay,
      name = "GLEAM_AMS",
      nbin = nbin,
      ctrl = ctrl,
      seed = seed
    )
    md <- out[[]]
    cols <- grep("^GLEAM_AMS", colnames(md), value = TRUE)
    if (length(cols) > length(genesets)) {
      cols <- utils::tail(cols, length(genesets))
    }
    mat <- t(as.matrix(md[, cols, drop = FALSE]))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(out)
    return(mat)
  }

  if (method == "UCell") {
    check_method_dependency(method, "UCell")
    expr <- as_expr_matrix(expr)
    dense <- as_dense_matrix(expr)
    max_rank <- method_params$ucell_max_rank %||% NULL
    u <- if (is.null(max_rank)) {
      UCell::ScoreSignatures_UCell(dense, features = genesets)
    } else {
      UCell::ScoreSignatures_UCell(dense, features = genesets, maxRank = max_rank)
    }
    mat <- t(as.matrix(u[, names(genesets), drop = FALSE]))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(mat)
  }

  if (method == "AUCell") {
    check_method_dependency(method, "AUCell")
    expr <- as_expr_matrix(expr)
    auc_max_rank <- method_params$auc_max_rank %||% 0.05
    auc_rank_abs <- if (is.numeric(auc_max_rank) && length(auc_max_rank) == 1L && auc_max_rank <= 1) {
      max(1L, ceiling(nrow(expr) * auc_max_rank))
    } else {
      as.integer(round(auc_max_rank))
    }
    rankings <- AUCell::AUCell_buildRankings(as.matrix(expr), plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(genesets, rankings, aucMaxRank = auc_rank_abs)
    mat <- as.matrix(AUCell::getAUC(auc))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(mat)
  }

  if (method == "ssGSEA") {
    check_method_dependency(method, "GSVA")
    expr <- as_expr_matrix(expr)
    alpha <- method_params$alpha %||% 0.25
    normalize <- method_params$normalize %||% TRUE
    gsva_ns <- asNamespace("GSVA")
    ssgsea_param_fn <- get0("ssgseaParam", envir = gsva_ns, mode = "function")

    mat <- if (!is.null(ssgsea_param_fn)) {
      param <- ssgsea_param_fn(exprData = as.matrix(expr), geneSets = genesets, alpha = alpha, normalize = normalize)
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(as.matrix(expr), genesets, method = "ssgsea", verbose = FALSE)
    }
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(as.matrix(mat))
  }

  if (method == "GSVA") {
    check_method_dependency(method, "GSVA")
    expr <- as_expr_matrix(expr)
    kcdf <- method_params$kcdf %||% "Gaussian"
    gsva_ns <- asNamespace("GSVA")
    gsva_param_fn <- get0("gsvaParam", envir = gsva_ns, mode = "function")

    mat <- if (!is.null(gsva_param_fn)) {
      param <- gsva_param_fn(exprData = as.matrix(expr), geneSets = genesets, kcdf = kcdf)
      GSVA::gsva(param, verbose = FALSE)
    } else {
      GSVA::gsva(as.matrix(expr), genesets, method = "gsva", kcdf = kcdf, verbose = FALSE)
    }
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(as.matrix(mat))
  }

  if (method == "singscore") {
    check_method_dependency(method, "singscore")
    expr <- as_expr_matrix(expr)
    out <- lapply(names(genesets), function(nm) {
      singscore::simpleScore(as.matrix(expr), upSet = genesets[[nm]])$TotalScore
    })
    mat <- do.call(rbind, out)
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(mat)
  }

  stop(sprintf("Optional scoring method '%s' is not implemented.", method), call. = FALSE)
}
