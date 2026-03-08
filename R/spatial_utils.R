#' Join score with spatial coordinates
#'
#' @param score `gleam_score` object.
#' @param coords Spatial coordinates data.frame or matrix with x/y.
#' @param meta Optional metadata to merge.
#'
#' @return Data.frame containing cell_id, coordinates, and metadata.
#' @export
join_score_spatial <- function(score, coords, meta = NULL) {
  check_score_object(score)
  coords <- as.data.frame(coords)
  check_required_columns(coords, c("x", "y"))

  if (is.null(rownames(coords))) {
    if (nrow(coords) != ncol(score$score)) {
      stop("If coordinates have no rownames, nrow(coords) must equal ncol(score).", call. = FALSE)
    }
    rownames(coords) <- colnames(score$score)
  }

  out <- data.frame(
    cell_id = colnames(score$score),
    x = coords[colnames(score$score), "x"],
    y = coords[colnames(score$score), "y"],
    stringsAsFactors = FALSE,
    row.names = colnames(score$score)
  )

  meta_use <- if (is.null(meta)) score$meta else as.data.frame(meta, stringsAsFactors = FALSE)
  if (!is.null(rownames(meta_use))) {
    meta_use <- meta_use[colnames(score$score), , drop = FALSE]
  }
  # Avoid duplicate columns (`cell_id`, `x`, `y`) when metadata already carries spatial fields.
  meta_use <- meta_use[, setdiff(colnames(meta_use), c("cell_id", "x", "y")), drop = FALSE]

  cbind(out, meta_use)
}

#' Test spatial region/domain pathway differences
#'
#' @param score `gleam_score` object.
#' @param region Region/domain metadata column or vector.
#' @param group Optional group for within-region contrasts.
#' @param sample Optional sample variable.
#' @param method Statistical method.
#' @param adjust_method Multiple testing adjustment.
#' @param level Region level (`region` or `sample_region`).
#'
#' @return `gleam_test` object.
#' @export
test_pathway_spatial <- function(
  score,
  region,
  group = NULL,
  sample = NULL,
  method = c("wilcox", "t"),
  adjust_method = "BH",
  level = c("region", "sample_region")
) {
  check_score_object(score)
  method <- match.arg(method)
  level <- match.arg(level)

  meta <- score$meta
  region_v <- resolve_meta_var(meta, region, "region")

  if (is.null(group)) {
    # Region vs region comparison.
    tbl <- test_pathway_celltype(
      score_mat = score$score,
      celltype = region_v,
      method = method,
      ref_celltype = NULL,
      adjust_method = adjust_method
    )
    tbl$comparison_type <- "region"
    tbl$level <- "region"
  } else {
    group_v <- resolve_meta_var(meta, group, "group")
    if (level == "region") {
      res <- lapply(unique(region_v), function(rg) {
        sub <- test_groups_within_celltype(
          score_mat = score$score,
          group = group_v,
          celltype = region_v,
          target_celltype = rg,
          level = "cell",
          method = method,
          adjust_method = adjust_method
        )
        sub$comparison_type <- "group_within_region"
        sub$level <- "region"
        sub$celltype <- rg
        sub
      })
      tbl <- do.call(rbind, res)
    } else {
      if (is.null(sample)) stop("`sample` is required for level = 'sample_region'.", call. = FALSE)
      sample_v <- resolve_meta_var(meta, sample, "sample")
      tbl <- test_pathway_sample_celltype(
        score_mat = score$score,
        group = group_v,
        sample = sample_v,
        celltype = region_v,
        method = method,
        adjust_method = adjust_method
      )
      tbl$comparison_type <- "group_within_region"
      tbl$level <- "sample_region"
    }
  }

  new_scpathway_test(
    table = tbl,
    level = level,
    method = method,
    comparison = list(type = "spatial", level = level),
    params = list(adjust_method = adjust_method)
  )
}
