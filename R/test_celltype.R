#' Celltype-level pathway test
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param celltype Celltype vector.
#' @param method Test method.
#' @param ref_celltype Reference celltype.
#' @param adjust_method P-value adjustment method.
#'
#' @return Differential pathway result data.frame.
#' @keywords internal
test_pathway_celltype <- function(
  score_mat,
  celltype,
  method = c("wilcox", "t"),
  ref_celltype = NULL,
  adjust_method = "BH"
) {
  method <- match.arg(method)
  test_pathway_cell(
    score_mat = score_mat,
    group = as.character(celltype),
    method = method,
    ref_group = ref_celltype,
    adjust_method = adjust_method,
    comparison_type = "celltype",
    celltype_label = NA_character_,
    level = "celltype"
  )
}

#' Sample+celltype-level pathway test
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param group Group vector.
#' @param sample Sample vector.
#' @param celltype Celltype vector.
#' @param method Test method.
#' @param ref_group Reference group.
#' @param adjust_method P-value adjustment method.
#'
#' @return Differential pathway result data.frame.
#' @keywords internal
test_pathway_sample_celltype <- function(
  score_mat,
  group,
  sample,
  celltype,
  method = c("wilcox", "t"),
  ref_group = NULL,
  adjust_method = "BH"
) {
  method <- match.arg(method)
  ct_levels <- unique(as.character(celltype))

  res <- lapply(ct_levels, function(ct) {
    idx <- as.character(celltype) == ct
    if (sum(idx) < 4L) {
      return(NULL)
    }
    if (length(unique(as.character(group)[idx])) < 2L) {
      return(NULL)
    }

    sub <- test_pathway_sample(
      score_mat = score_mat[, idx, drop = FALSE],
      group = as.character(group)[idx],
      sample = as.character(sample)[idx],
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method
    )
    sub$celltype <- ct
    sub$comparison_type <- "group_within_celltype"
    sub$level <- "sample_celltype"
    sub
  })

  res <- res[!vapply(res, is.null, logical(1))]
  if (length(res) == 0L) {
    return(data.frame(
      pathway = character(0),
      comparison_type = character(0),
      group1 = character(0),
      group2 = character(0),
      celltype = character(0),
      level = character(0),
      effect_size = numeric(0),
      median_group1 = numeric(0),
      median_group2 = numeric(0),
      diff_median = numeric(0),
      p_value = numeric(0),
      p_adj = numeric(0),
      n_group1 = integer(0),
      n_group2 = integer(0),
      mean_group1 = numeric(0),
      mean_group2 = numeric(0),
      direction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}
