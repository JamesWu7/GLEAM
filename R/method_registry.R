#' List available scoring methods
#'
#' @param include_optional Whether to include optional integration methods.
#'
#' @return A data.frame describing methods and dependencies.
#' @export
list_scoring_methods <- function(include_optional = TRUE) {
  reg <- get_scoring_registry()
  if (!include_optional) {
    reg <- reg[reg$category == "native", , drop = FALSE]
  }
  rownames(reg) <- NULL
  reg
}

#' @keywords internal
get_scoring_registry <- function() {
  data.frame(
    method_name = c(
      "rank", "auc", "zmean", "mean", "scaled_mean", "robust", "singscore_like", "ssgsea_like", "ensemble",
      "addmodulescore", "ucell", "aucell", "gsva", "singscore"
    ),
    category = c(rep("native", 9), rep("optional", 5)),
    dependency = c(
      rep(NA_character_, 9),
      "Seurat", "UCell", "AUCell", "GSVA", "singscore"
    ),
    description = c(
      "Rank-based normalized score",
      "Top-rank AUC-like score",
      "Mean standardized expression",
      "Mean expression",
      "Z-scored mean expression",
      "Robust median/MAD normalized score",
      "Singscore-inspired rank score",
      "ssGSEA-like rank enrichment approximation",
      "Combine multiple methods",
      "Seurat AddModuleScore integration",
      "UCell integration",
      "AUCell integration",
      "GSVA/ssGSEA integration",
      "singscore integration"
    ),
    input_assumption = c(
      rep("gene_by_cell matrix", 9),
      "Seurat object", "gene_by_cell matrix", "gene_by_cell matrix", "gene_by_cell matrix", "gene_by_cell matrix"
    ),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
check_scoring_method <- function(method) {
  reg <- get_scoring_registry()
  if (!method %in% reg$method_name) {
    stop(sprintf("Unsupported scoring method '%s'. Use list_scoring_methods().", method), call. = FALSE)
  }
  row <- reg[reg$method_name == method, , drop = FALSE]
  if (row$category[[1]] == "optional") {
    check_method_dependency(method = method, pkg = row$dependency[[1]])
  }
  invisible(TRUE)
}
