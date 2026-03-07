#' Volcano plot for differential pathways
#'
#' @param x `scpathway_test` object or result data.frame.
#' @param p_col P-value column name.
#' @param effect_col Effect size column name.
#' @param sig_thresh Significance threshold for adjusted p-value.
#'
#' @return A `ggplot` object.
#' @export
plot_volcano <- function(x, p_col = "p_adj", effect_col = "effect_size", sig_thresh = 0.05) {
  tbl <- if (inherits(x, "scpathway_test")) x$table else x
  if (!is.data.frame(tbl)) {
    stop("`x` must be a scpathway_test object or a data.frame.", call. = FALSE)
  }

  if (!all(c(p_col, effect_col) %in% colnames(tbl))) {
    stop("Missing required columns for volcano plot.", call. = FALSE)
  }

  pval <- tbl[[p_col]]
  pval[is.na(pval) | pval <= 0] <- 1
  y <- -log10(pval)
  sig <- tbl[[p_col]] < sig_thresh

  df <- data.frame(
    effect = tbl[[effect_col]],
    neglog10 = y,
    significant = ifelse(sig, "yes", "no"),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = effect, y = neglog10, color = significant)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::geom_hline(yintercept = -log10(sig_thresh), linetype = 2, color = "grey40") +
    ggplot2::labs(x = "Effect size", y = "-log10(adj p)", title = "Differential pathway volcano") +
    .theme_scpathway()
}
