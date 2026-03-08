#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
stopf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

#' @keywords internal
warnf <- function(fmt, ...) {
  warning(sprintf(fmt, ...), call. = FALSE)
}

#' @keywords internal
check_required_columns <- function(meta, cols) {
  miss <- setdiff(cols, colnames(meta))
  if (length(miss) > 0) {
    stopf("Missing metadata columns: %s", paste(miss, collapse = ", "))
  }
}

#' @keywords internal
resolve_meta_var <- function(meta, x, arg_name) {
  if (length(x) == 1L && is.character(x) && x %in% colnames(meta)) {
    return(meta[[x]])
  }
  if (length(x) == nrow(meta)) {
    return(x)
  }
  stopf("`%s` must be a metadata column name or a vector of length ncol(score).", arg_name)
}

#' @keywords internal
check_score_object <- function(score) {
  if (!inherits(score, "gleam_score") && !inherits(score, "scpathway_score")) {
    stop("`score` must be an object returned by score_pathway().", call. = FALSE)
  }
}

#' @keywords internal
require_optional_package <- function(pkg, feature = NULL) {
  if (identical(pkg, "monocle3")) {
    return(.assert_monocle3(feature = feature %||% "trajectory analysis"))
  }
  .assert_pkg(pkg = pkg, feature = feature)
}

#' @keywords internal
check_method_dependency <- function(method, pkg = NULL) {
  if (!is.null(pkg)) {
    require_optional_package(pkg, feature = sprintf("method '%s'", method))
  }
  invisible(TRUE)
}

#' @keywords internal
check_plot_dependency <- function(plot_name, pkg) {
  require_optional_package(pkg, feature = sprintf("plot '%s'", plot_name))
  invisible(TRUE)
}

#' @keywords internal
match_arg_chr <- function(arg, choices, arg_name) {
  arg <- match.arg(arg, choices)
  if (length(arg) != 1L) {
    stopf("`%s` must resolve to a single value.", arg_name)
  }
  arg
}
