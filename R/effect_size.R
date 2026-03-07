#' Median difference effect size
#'
#' @param x Numeric score vector.
#' @param g Group vector with two groups.
#'
#' @return Numeric effect size.
#' @keywords internal
effect_median_diff <- function(x, g) {
  g <- as.character(g)
  lev <- unique(g[!is.na(g)])
  if (length(lev) != 2L) return(NA_real_)
  stats::median(x[g == lev[1]], na.rm = TRUE) - stats::median(x[g == lev[2]], na.rm = TRUE)
}

#' Rank-biserial effect size placeholder
#'
#' @param x Numeric vector.
#' @param g Group vector.
#'
#' @return Numeric or NA.
#' @keywords internal
effect_rank_biserial <- function(x, g) {
  NA_real_
}

#' Standardized mean difference placeholder
#'
#' @param x Numeric vector.
#' @param g Group vector.
#'
#' @return Numeric or NA.
#' @keywords internal
effect_standardized_diff <- function(x, g) {
  g <- as.character(g)
  lev <- unique(g[!is.na(g)])
  if (length(lev) != 2L) return(NA_real_)
  x1 <- x[g == lev[1]]
  x2 <- x[g == lev[2]]
  s <- stats::sd(c(x1, x2), na.rm = TRUE)
  if (is.na(s) || s == 0) return(NA_real_)
  (mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)) / s
}
