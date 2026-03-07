#' Violin plot of pathway scores
#'
#' @param score `scpathway_score` object.
#' @param pathway Pathway name.
#' @param group Group variable (metadata column name or vector).
#' @param celltype Optional celltype filter (single label or vector).
#' @param trim Passed to `geom_violin()`.
#'
#' @return A `ggplot` object.
#' @export
plot_violin <- function(score, pathway, group, celltype = NULL, trim = TRUE) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) {
    stop("`pathway` not found in score matrix.", call. = FALSE)
  }

  meta <- score$meta
  g <- resolve_meta_var(meta, group, "group")
  y <- as.numeric(score$score[pathway, ])

  keep <- rep(TRUE, length(y))
  if (!is.null(celltype)) {
    ct <- if (length(celltype) == 1L && is.character(celltype) && celltype %in% colnames(meta)) {
      as.character(meta[[celltype]])
    } else if (length(celltype) == 1L) {
      as.character(meta$celltype)
    } else {
      as.character(celltype)
    }
    if (length(ct) == 1L) {
      keep <- meta$celltype %in% ct
    }
  }

  df <- data.frame(group = as.factor(g[keep]), score = y[keep], stringsAsFactors = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = group, y = score, fill = group)) +
    ggplot2::geom_violin(trim = trim, alpha = 0.7, color = NA) +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8) +
    ggplot2::labs(x = "Group", y = pathway, title = paste("Pathway score:", pathway)) +
    .theme_scpathway()
}
