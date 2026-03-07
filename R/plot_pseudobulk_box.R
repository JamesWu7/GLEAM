#' Pseudobulk boxplot for pathway scores
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param group Group variable.
#' @param sample Sample variable.
#' @param celltype Optional celltype variable.
#'
#' @return A `ggplot` object.
#' @export
plot_pseudobulk_box <- function(score, pathway, group, sample, celltype = NULL) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)

  g <- resolve_meta_var(score$meta, group, "group")
  s <- resolve_meta_var(score$meta, sample, "sample")
  ctv <- if (!is.null(celltype)) resolve_meta_var(score$meta, celltype, "celltype") else rep("all", ncol(score$score))

  df <- data.frame(
    sample = as.character(s),
    group = as.character(g),
    celltype = as.character(ctv),
    value = as.numeric(score$score[pathway, ]),
    stringsAsFactors = FALSE
  )

  pb <- stats::aggregate(value ~ sample + group + celltype, data = df, FUN = mean)
  pb$group_celltype <- interaction(pb$group, pb$celltype, drop = TRUE)

  ggplot2::ggplot(pb, ggplot2::aes(x = ggplot2::.data$group_celltype, y = ggplot2::.data$value, fill = ggplot2::.data$group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    ggplot2::geom_jitter(width = 0.12, size = 1.6, alpha = 0.9) +
    ggplot2::labs(title = paste("Pseudobulk:", pathway), x = "Group-Celltype", y = "Mean pathway score") +
    .theme_gleam()
}
