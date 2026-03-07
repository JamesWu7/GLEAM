#' @keywords internal
run_group_test <- function(x, g, method = c("wilcox", "t")) {
  method <- match.arg(method)
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- as.character(g[keep])

  lev <- unique(g)
  if (length(lev) != 2L) {
    return(list(p_value = NA_real_, n1 = NA_integer_, n2 = NA_integer_))
  }

  x1 <- x[g == lev[1]]
  x2 <- x[g == lev[2]]

  if (length(x1) < 2L || length(x2) < 2L) {
    return(list(p_value = NA_real_, n1 = length(x1), n2 = length(x2)))
  }

  p <- tryCatch(
    if (method == "wilcox") {
      suppressWarnings(stats::wilcox.test(x1, x2)$p.value)
    } else {
      stats::t.test(x1, x2)$p.value
    },
    error = function(e) NA_real_
  )

  list(p_value = as.numeric(p), n1 = length(x1), n2 = length(x2))
}

#' @keywords internal
compute_effects <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]
  g <- as.character(g[keep])
  lev <- unique(g)
  if (length(lev) != 2L) {
    return(list(
      group1 = NA_character_, group2 = NA_character_, median_group1 = NA_real_,
      median_group2 = NA_real_, diff_median = NA_real_, effect_size = NA_real_,
      mean_group1 = NA_real_, mean_group2 = NA_real_
    ))
  }

  x1 <- x[g == lev[1]]
  x2 <- x[g == lev[2]]
  med1 <- stats::median(x1, na.rm = TRUE)
  med2 <- stats::median(x2, na.rm = TRUE)

  list(
    group1 = lev[1],
    group2 = lev[2],
    median_group1 = med1,
    median_group2 = med2,
    diff_median = med1 - med2,
    effect_size = med1 - med2,
    mean_group1 = mean(x1, na.rm = TRUE),
    mean_group2 = mean(x2, na.rm = TRUE)
  )
}

#' @keywords internal
pick_two_groups <- function(g, ref_group = NULL) {
  g <- as.character(g)
  lev <- sort(unique(g[!is.na(g)]))
  if (length(lev) < 2L) {
    stop("Need at least two groups for comparison.", call. = FALSE)
  }
  if (!is.null(ref_group)) {
    if (!ref_group %in% lev) {
      stop("`ref_group` is not present in groups.", call. = FALSE)
    }
    other <- setdiff(lev, ref_group)
    if (length(other) < 1L) {
      stop("No comparison group available after applying `ref_group`.", call. = FALSE)
    }
    return(c(ref_group, other[1]))
  }
  lev[1:2]
}

#' @keywords internal
p_adjust_safe <- function(p, method = "BH") {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], method = method)
  out
}
