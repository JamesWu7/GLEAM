#' @keywords internal
score_optional_wrapper <- function(method, object = NULL, expr = NULL, genesets, seurat = FALSE, assay = NULL, slot = "data") {
  method <- tolower(method)

  if (method == "addmodulescore") {
    check_method_dependency(method, "Seurat")
    if (!seurat || is.null(object)) {
      stop("`addmodulescore` requires `seurat = TRUE` and a Seurat object.", call. = FALSE)
    }
    nm <- paste0("GLEAM_AMS_", seq_along(genesets))
    out <- Seurat::AddModuleScore(object = object, features = unname(genesets), assay = assay, name = "GLEAM_AMS")
    md <- out[[]]
    cols <- grep("^GLEAM_AMS", colnames(md), value = TRUE)
    mat <- t(as.matrix(md[, cols, drop = FALSE]))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(out)
    return(mat)
  }

  if (method == "ucell") {
    check_method_dependency(method, "UCell")
    expr <- as_expr_matrix(expr)
    dense <- as_dense_matrix(expr)
    u <- UCell::ScoreSignatures_UCell(dense, features = genesets)
    mat <- t(as.matrix(u[, names(genesets), drop = FALSE]))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(mat)
  }

  if (method == "aucell") {
    check_method_dependency(method, "AUCell")
    expr <- as_expr_matrix(expr)
    rankings <- AUCell::AUCell_buildRankings(as.matrix(expr), plotStats = FALSE, verbose = FALSE)
    auc <- AUCell::AUCell_calcAUC(genesets, rankings)
    mat <- as.matrix(AUCell::getAUC(auc))
    rownames(mat) <- names(genesets)
    colnames(mat) <- colnames(expr)
    return(mat)
  }

  if (method == "gsva") {
    check_method_dependency(method, "GSVA")
    expr <- as_expr_matrix(expr)
    param <- GSVA::ssgseaParam(expr = as.matrix(expr), geneSets = genesets)
    mat <- GSVA::gsva(param)
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
