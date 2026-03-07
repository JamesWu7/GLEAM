#' Test pathway differences
#'
#' @param score `scpathway_score` object.
#' @param group Group variable (column name or vector).
#' @param sample Sample variable (column name or vector).
#' @param celltype Celltype variable (column name or vector).
#' @param level Comparison level.
#' @param method Statistical method.
#' @param ref_group Reference group.
#' @param paired Placeholder, currently ignored.
#' @param covariates Placeholder, currently ignored.
#' @param adjust_method Multiple testing adjustment method.
#' @param min_cells Minimum cells per group (reserved for future checks).
#' @param min_samples Minimum samples per group (reserved for future checks).
#' @param verbose Print messages.
#'
#' @return An object of class `scpathway_test`.
#' @export
test_pathway <- function(
  score,
  group,
  sample = NULL,
  celltype = NULL,
  level = c("cell", "sample", "celltype", "sample_celltype"),
  method = c("wilcox", "t", "lm"),
  ref_group = NULL,
  paired = FALSE,
  covariates = NULL,
  adjust_method = "BH",
  min_cells = 10,
  min_samples = 2,
  verbose = TRUE
) {
  check_score_object(score)
  level <- match.arg(level)
  method <- match.arg(method)
  if (method == "lm") {
    warning("lm method is reserved/experimental in v1.", call. = FALSE)
    method <- "wilcox"
  }

  meta <- score$meta
  group_v <- resolve_meta_var(meta, group, "group")

  if (level == "cell") {
    if (verbose) warning("Cell-level testing is exploratory. Prefer sample-level for formal inference.", call. = FALSE)
    tbl <- test_pathway_cell(
      score_mat = score$score,
      group = group_v,
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method,
      level = "cell"
    )
  } else if (level == "sample") {
    if (is.null(sample)) stop("`sample` is required for level = 'sample'.", call. = FALSE)
    sample_v <- resolve_meta_var(meta, sample, "sample")
    tbl <- test_pathway_sample(
      score_mat = score$score,
      group = group_v,
      sample = sample_v,
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method
    )
  } else if (level == "celltype") {
    if (is.null(celltype)) stop("`celltype` is required for level = 'celltype'.", call. = FALSE)
    celltype_v <- resolve_meta_var(meta, celltype, "celltype")
    tbl <- test_pathway_celltype(
      score_mat = score$score,
      celltype = celltype_v,
      method = method,
      ref_celltype = ref_group,
      adjust_method = adjust_method
    )
  } else {
    if (is.null(sample) || is.null(celltype)) {
      stop("`sample` and `celltype` are required for level = 'sample_celltype'.", call. = FALSE)
    }
    sample_v <- resolve_meta_var(meta, sample, "sample")
    celltype_v <- resolve_meta_var(meta, celltype, "celltype")
    tbl <- test_pathway_sample_celltype(
      score_mat = score$score,
      group = group_v,
      sample = sample_v,
      celltype = celltype_v,
      method = method,
      ref_group = ref_group,
      adjust_method = adjust_method
    )
  }

  new_scpathway_test(
    table = tbl,
    level = level,
    method = method,
    comparison = list(level = level),
    params = list(
      paired = paired,
      covariates = covariates,
      adjust_method = adjust_method,
      min_cells = min_cells,
      min_samples = min_samples
    )
  )
}

#' Compare pathways across cell types
#'
#' @param score `scpathway_score` object.
#' @param celltype Celltype variable name or vector.
#' @param group Optional group variable (reserved for future stratified mode).
#' @param method Statistical method.
#' @param ref_celltype Reference celltype.
#' @param adjust_method Multiple testing adjustment method.
#' @param verbose Print messages.
#'
#' @return An object of class `scpathway_test`.
#' @export
compare_celltypes <- function(
  score,
  celltype,
  group = NULL,
  method = c("wilcox", "t"),
  ref_celltype = NULL,
  adjust_method = "BH",
  verbose = TRUE
) {
  method <- match.arg(method)
  if (!is.null(group) && verbose) {
    message("`group` argument is currently reserved for future stratified celltype comparisons.")
  }
  test_pathway(
    score = score,
    group = if (is.null(group)) rep("all", ncol(score$score)) else group,
    celltype = celltype,
    level = "celltype",
    method = method,
    ref_group = ref_celltype,
    adjust_method = adjust_method,
    verbose = verbose
  )
}

#' Compare groups within a fixed celltype
#'
#' @param score `scpathway_score` object.
#' @param group Group variable name or vector.
#' @param celltype Celltype variable name or vector.
#' @param target_celltype Target celltype label.
#' @param sample Optional sample variable for sample-level testing.
#' @param level Comparison level: `cell` or `sample`.
#' @param method Statistical method.
#' @param ref_group Optional reference group.
#' @param adjust_method Multiple testing adjustment method.
#' @param verbose Print messages.
#'
#' @return An object of class `scpathway_test`.
#' @export
compare_groups_within_celltype <- function(
  score,
  group,
  celltype,
  target_celltype,
  sample = NULL,
  level = c("cell", "sample"),
  method = c("wilcox", "t", "lm"),
  ref_group = NULL,
  adjust_method = "BH",
  verbose = TRUE
) {
  check_score_object(score)
  level <- match.arg(level)
  method <- match.arg(method)
  if (method == "lm") {
    warning("lm method is reserved/experimental in v1, falling back to wilcox.", call. = FALSE)
    method <- "wilcox"
  }

  meta <- score$meta
  g <- resolve_meta_var(meta, group, "group")
  ct <- resolve_meta_var(meta, celltype, "celltype")
  smp <- if (!is.null(sample)) resolve_meta_var(meta, sample, "sample") else NULL

  tbl <- test_groups_within_celltype(
    score_mat = score$score,
    group = g,
    celltype = ct,
    target_celltype = target_celltype,
    sample = smp,
    level = level,
    method = if (method == "t") "t" else "wilcox",
    ref_group = ref_group,
    adjust_method = adjust_method
  )

  new_scpathway_test(
    table = tbl,
    level = level,
    method = method,
    comparison = list(target_celltype = target_celltype),
    params = list(adjust_method = adjust_method)
  )
}
