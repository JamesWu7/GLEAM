#' Sample-level pathway test
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param group Group vector.
#' @param sample Sample vector.
#' @param method Test method.
#' @param ref_group Reference group.
#' @param adjust_method P-value adjustment method.
#'
#' @return Differential pathway result data.frame.
#' @keywords internal
test_pathway_sample <- function(
  score_mat,
  group,
  sample,
  method = c("wilcox", "t"),
  ref_group = NULL,
  adjust_method = "BH"
) {
  method <- match.arg(method)
  df <- data.frame(
    cell = colnames(score_mat),
    sample = as.character(sample),
    group = as.character(group),
    stringsAsFactors = FALSE
  )

  samp_levels <- unique(df$sample)
  if (length(samp_levels) < 2L) {
    stop("Need at least two samples for sample-level testing.", call. = FALSE)
  }

  sample_group <- tapply(df$group, df$sample, function(x) unique(x)[1])
  sample_score <- sapply(samp_levels, function(s) {
    idx <- df$sample == s
    rowMeans(score_mat[, idx, drop = FALSE], na.rm = TRUE)
  })

  group_s <- as.character(sample_group[colnames(sample_score)])
  test_pathway_cell(
    score_mat = sample_score,
    group = group_s,
    method = method,
    ref_group = ref_group,
    adjust_method = adjust_method,
    comparison_type = "group",
    celltype_label = NA_character_,
    level = "sample"
  )
}
