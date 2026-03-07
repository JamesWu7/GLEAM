#' Plot pathway score over pseudotime
#'
#' @param score `gleam_score` object.
#' @param pathway Pathway name.
#' @param pseudotime Pseudotime source.
#' @param lineage Optional lineage source for coloring.
#' @param smooth Add smoothing line.
#'
#' @return A `ggplot` object.
#' @export
plot_pseudotime_score <- function(score, pathway, pseudotime = NULL, lineage = NULL, smooth = TRUE) {
  check_score_object(score)
  if (!pathway %in% rownames(score$score)) stop("`pathway` not found in score matrix.", call. = FALSE)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  df <- data.frame(
    pseudotime = pt,
    score = as.numeric(score$score[pathway, ]),
    lineage = as.factor(ln),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = ggplot2::.data$pseudotime, y = ggplot2::.data$score, color = ggplot2::.data$lineage)) +
    ggplot2::geom_point(alpha = 0.55, size = 1.1) +
    ggplot2::labs(title = paste("Pseudotime score:", pathway), x = "Pseudotime", y = "Pathway score") +
    .theme_gleam()

  if (smooth) {
    p <- p + ggplot2::geom_smooth(se = FALSE, method = "loess", formula = y ~ x)
  }
  p
}
