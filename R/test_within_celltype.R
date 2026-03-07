#' Within-celltype group comparison
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param group Group vector.
#' @param celltype Celltype vector.
#' @param target_celltype Target celltype.
#' @param sample Optional sample vector.
#' @param level Comparison level (`cell` or `sample`).
#' @param method Test method.
#' @param ref_group Reference group.
#' @param adjust_method P-value adjustment method.
#'
#' @return Differential pathway result data.frame.
#' @keywords internal
test_groups_within_celltype <- function(
  score_mat,
  group,
  celltype,
  target_celltype,
  sample = NULL,
  level = c("cell", "sample"),
  method = c("wilcox", "t"),
  ref_group = NULL,
  adjust_method = "BH"
) {
  level <- match.arg(level)
  method <- match.arg(method)

  idx <- as.character(celltype) == target_celltype
  if (!any(idx)) {
    stop("No cells matched `target_celltype`.", call. = FALSE)
  }

  mat <- score_mat[, idx, drop = FALSE]
  g <- as.character(group)[idx]

  if (level == "cell") {
    out <- test_pathway_cell(
      score_mat = mat,
      group = g,
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method,
      comparison_type = "group_within_celltype",
      celltype_label = target_celltype,
      level = "cell"
    )
  } else {
    if (is.null(sample)) {
      stop("`sample` must be provided for level = 'sample'.", call. = FALSE)
    }
    out <- test_pathway_sample(
      score_mat = mat,
      group = g,
      sample = as.character(sample)[idx],
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method
    )
    out$comparison_type <- "group_within_celltype"
    out$celltype <- target_celltype
    out$level <- "sample"
  }

  out
}
