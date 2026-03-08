#' Test pathway differences
#'
#' @param score `gleam_score` object.
#' @param group Group variable (column name or vector).
#' @param sample Sample variable (column name or vector).
#' @param celltype Celltype variable (column name or vector).
#' @param target_celltype Target celltype label for within-celltype group comparisons.
#' @param region Spatial region/domain variable (column name or vector).
#' @param pseudotime Pseudotime source for trajectory mode.
#' @param lineage Lineage source for trajectory mode.
#' @param level Comparison level.
#' @param method Statistical method.
#' @param ref_group Reference group.
#' @param paired Placeholder, currently ignored.
#' @param covariates Placeholder, currently ignored.
#' @param adjust_method Multiple testing adjustment method.
#' @param min_cells Minimum cells per group.
#' @param min_samples Minimum samples per group.
#' @param aggregation Aggregation summary used by pseudobulk levels.
#' @param threshold Threshold used when `aggregation = 'fraction'`.
#' @param backend Trajectory backend used when `level = 'trajectory'`.
#' @param verbose Print messages.
#'
#' @return An object of class `gleam_test`.
#' @keywords internal
#' @noRd
test_pathway <- function(
  score,
  group = NULL,
  sample = NULL,
  celltype = NULL,
  target_celltype = NULL,
  region = NULL,
  pseudotime = NULL,
  lineage = NULL,
  level = c("cell", "sample", "celltype", "sample_celltype", "pseudobulk", "region", "sample_region", "trajectory"),
  method = c("wilcox", "t", "lm", "limma", "edgeR", "DESeq2", "spearman", "tradeSeq"),
  ref_group = NULL,
  paired = FALSE,
  covariates = NULL,
  adjust_method = "BH",
  min_cells = 10,
  min_samples = 2,
  aggregation = c("mean", "median", "fraction", "sum"),
  threshold = 0,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot"),
  verbose = TRUE
) {
  check_score_object(score)
  level <- match.arg(level)
  method <- match.arg(method)
  aggregation <- match.arg(aggregation)
  backend <- match.arg(backend)

  meta <- score$meta
  if (is.null(group) && level %in% c("cell", "sample", "sample_celltype", "pseudobulk", "sample_region")) {
    stop("`group` is required for selected comparison level.", call. = FALSE)
  }
  group_v <- if (!is.null(group)) resolve_meta_var(meta, group, "group") else NULL

  if (!is.null(target_celltype)) {
    if (is.null(group) || is.null(celltype)) {
      stop("`group` and `celltype` are required when `target_celltype` is provided.", call. = FALSE)
    }
    use_level <- if (level %in% c("cell", "sample")) level else if (!is.null(sample)) "sample" else "cell"
    smp <- if (!is.null(sample)) resolve_meta_var(meta, sample, "sample") else NULL
    tbl <- test_groups_within_celltype(
      score_mat = score$score,
      group = group_v,
      celltype = resolve_meta_var(meta, celltype, "celltype"),
      target_celltype = target_celltype,
      sample = smp,
      level = use_level,
      method = if (method == "t") "t" else "wilcox",
      ref_group = ref_group,
      adjust_method = adjust_method
    )
    return(new_scpathway_test(
      table = tbl,
      level = use_level,
      method = method,
      comparison = list(target_celltype = target_celltype),
      params = list(
        paired = paired,
        covariates = covariates,
        adjust_method = adjust_method,
        min_cells = min_cells,
        min_samples = min_samples,
        aggregation = aggregation,
        threshold = threshold,
        backend = backend
      )
    ))
  }

  if (level == "trajectory") {
    tst <- test_pathway_trajectory(
      score = score,
      pathway = NULL,
      pseudotime = pseudotime,
      lineage = lineage,
      method = if (method %in% c("spearman", "lm", "tradeSeq")) method else "spearman",
      adjust_method = adjust_method,
      backend = backend,
      verbose = verbose
    )
    return(tst)
  }

  if (level == "cell") {
    if (verbose) warning("Cell-level testing is exploratory. Prefer sample/pseudobulk-level for formal inference.", call. = FALSE)
    tbl <- test_pathway_cell(
      score_mat = score$score,
      group = group_v,
      method = if (method == "t") "t" else "wilcox",
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
      method = if (method == "t") "t" else "wilcox",
      ref_group = ref_group,
      adjust_method = adjust_method
    )
  } else if (level == "celltype") {
    if (is.null(celltype)) stop("`celltype` is required for level = 'celltype'.", call. = FALSE)
    celltype_v <- resolve_meta_var(meta, celltype, "celltype")
    tbl <- test_pathway_celltype(
      score_mat = score$score,
      celltype = celltype_v,
      method = if (method == "t") "t" else "wilcox",
      ref_celltype = ref_group,
      adjust_method = adjust_method
    )
  } else if (level == "sample_celltype") {
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
      method = if (method == "t") "t" else "wilcox",
      ref_group = ref_group,
      adjust_method = adjust_method
    )
  } else if (level == "pseudobulk") {
    if (is.null(sample)) stop("`sample` is required for level = 'pseudobulk'.", call. = FALSE)
    sample_v <- resolve_meta_var(meta, sample, "sample")
    celltype_v <- if (!is.null(celltype)) resolve_meta_var(meta, celltype, "celltype") else NULL
    region_v <- if (!is.null(region)) resolve_meta_var(meta, region, "region") else NULL

    tbl <- test_pathway_pseudobulk(
      score = score,
      group = group_v,
      sample = sample_v,
      celltype = celltype_v,
      region = region_v,
      method = method,
      adjust_method = adjust_method,
      aggregation = aggregation,
      threshold = threshold,
      ref_group = ref_group,
      level = "pseudobulk"
    )
  } else if (level == "region") {
    if (is.null(region)) stop("`region` is required for level = 'region'.", call. = FALSE)
    region_v <- resolve_meta_var(meta, region, "region")
    tbl <- test_pathway_spatial(
      score = score,
      region = region_v,
      group = group_v,
      sample = sample,
      method = if (method == "t") "t" else "wilcox",
      adjust_method = adjust_method,
      level = "region"
    )$table
  } else { # sample_region
    if (is.null(sample) || is.null(region)) stop("`sample` and `region` are required for level = 'sample_region'.", call. = FALSE)
    sample_v <- resolve_meta_var(meta, sample, "sample")
    region_v <- resolve_meta_var(meta, region, "region")
    tbl <- test_pathway_pseudobulk(
      score = score,
      group = group_v,
      sample = sample_v,
      region = region_v,
      method = method,
      adjust_method = adjust_method,
      aggregation = aggregation,
      threshold = threshold,
      ref_group = ref_group,
      level = "sample_region"
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
      min_samples = min_samples,
      aggregation = aggregation,
      threshold = threshold,
      backend = backend
    )
  )
}

#' Compare pathways across cell types
#'
#' @param score `gleam_score` object.
#' @param celltype Celltype variable name or vector.
#' @param group Optional group variable (reserved for future stratified mode).
#' @param method Statistical method.
#' @param ref_celltype Reference celltype.
#' @param adjust_method Multiple testing adjustment method.
#' @param verbose Print messages.
#'
#' @return An object of class `gleam_test`.
#' @keywords internal
#' @noRd
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
#' @param score `gleam_score` object.
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
#' @return An object of class `gleam_test`.
#' @keywords internal
#' @noRd
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
    warning("lm method is reserved/experimental in v0.2, falling back to wilcox.", call. = FALSE)
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
