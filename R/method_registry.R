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
canonicalize_scoring_method <- function(method) {
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("`method` must be a non-empty character scalar.", call. = FALSE)
  }

  key <- tolower(trimws(as.character(method)))
  alias_map <- c(
    rank = "rank",
    mean = "mean",
    zscore = "zscore",
    scaled_mean = "scaled_mean",
    robust_mean = "robust_mean",
    ensemble = "ensemble",
    addmodulescore = "AddModuleScore",
    ucell = "UCell",
    aucell = "AUCell",
    ssgsea = "ssGSEA",
    gsva = "GSVA",
    singscore = "singscore"
  )

  out <- alias_map[[key]]
  if (is.null(out)) {
    stop(sprintf("Unsupported scoring method '%s'. Use list_scoring_methods().", method), call. = FALSE)
  }
  out
}

#' @keywords internal
canonicalize_scoring_methods <- function(methods) {
  methods <- vapply(methods, canonicalize_scoring_method, character(1))
  unique(unname(methods))
}

#' @keywords internal
get_scoring_registry <- function() {
  data.frame(
    method_name = c(
      "rank",
      "mean",
      "zscore",
      "scaled_mean",
      "robust_mean",
      "ensemble",
      "AddModuleScore",
      "UCell",
      "AUCell",
      "ssGSEA",
      "GSVA",
      "singscore"
    ),
    category = c(rep("native", 6), rep("optional", 6)),
    dependency = c(
      rep(NA_character_, 6),
      "Seurat", "UCell", "AUCell", "GSVA", "GSVA", "singscore"
    ),
    description = c(
      "Rank-based aggregate score (singscore-style concept)",
      "Average expression across signature genes",
      "Per-gene z-score followed by pathway aggregation",
      "Scaled mean score (equivalent to zscore aggregation)",
      "Robust aggregation using median/MAD scaling",
      "Ensemble of multiple methods after scale harmonization",
      "Seurat AddModuleScore integration",
      "UCell integration",
      "AUCell integration",
      "ssGSEA via GSVA package",
      "GSVA via GSVA package",
      "singscore integration"
    ),
    input_assumption = c(
      rep("gene_by_cell matrix", 6),
      "Seurat object", "gene_by_cell matrix", "gene_by_cell matrix", "gene_by_cell matrix", "gene_by_cell matrix", "gene_by_cell matrix"
    ),
    aliases = c(
      "",
      "",
      "",
      "",
      "",
      "",
      "addmodulescore",
      "ucell",
      "aucell",
      "",
      "gsva",
      ""
    ),
    is_legacy = rep(FALSE, 12),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
check_scoring_method <- function(method) {
  method <- canonicalize_scoring_method(method)
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
