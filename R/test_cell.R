#' Cell-level pathway test (exploratory)
#'
#' @param score_mat Pathway-by-cell matrix.
#' @param group Group vector.
#' @param method Test method.
#' @param ref_group Reference group.
#' @param adjust_method P-value adjustment method.
#' @param comparison_type Comparison type label.
#' @param celltype_label Optional celltype label.
#' @param level Result level label.
#'
#' @return Differential pathway result data.frame.
#' @keywords internal
test_pathway_cell <- function(
  score_mat,
  group,
  method = c("wilcox", "t"),
  ref_group = NULL,
  adjust_method = "BH",
  comparison_type = "group",
  celltype_label = NA_character_,
  level = "cell"
) {
  method <- match.arg(method)
  grp_pair <- pick_two_groups(group, ref_group = ref_group)
  keep <- as.character(group) %in% grp_pair
  g <- as.character(group)[keep]
  mat <- score_mat[, keep, drop = FALSE]

  rows <- lapply(seq_len(nrow(mat)), function(i) {
    x <- as.numeric(mat[i, ])
    eff <- compute_effects(x, g)
    tst <- run_group_test(x, g, method = method)

    data.frame(
      pathway = rownames(mat)[i],
      comparison_type = comparison_type,
      group1 = eff$group1,
      group2 = eff$group2,
      celltype = celltype_label,
      level = level,
      effect_size = eff$effect_size,
      median_group1 = eff$median_group1,
      median_group2 = eff$median_group2,
      diff_median = eff$diff_median,
      p_value = tst$p_value,
      p_adj = NA_real_,
      n_group1 = tst$n1,
      n_group2 = tst$n2,
      mean_group1 = eff$mean_group1,
      mean_group2 = eff$mean_group2,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$p_adj <- p_adjust_safe(out$p_value, method = adjust_method)
  out$direction <- ifelse(is.na(out$effect_size), NA_character_, ifelse(out$effect_size >= 0, "up", "down"))
  out
}
