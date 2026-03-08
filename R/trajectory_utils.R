#' Extract pseudotime values
#'
#' @param score `gleam_score` object.
#' @param pseudotime Metadata column name, numeric vector, or NULL.
#' @param object Optional trajectory object.
#'
#' @return Numeric vector.
#' @export
extract_pseudotime <- function(score, pseudotime = NULL, object = NULL) {
  check_score_object(score)
  meta <- score$meta

  if (!is.null(pseudotime) && is.null(object) && (inherits(pseudotime, "SlingshotDataSet") || inherits(pseudotime, "cell_data_set") || inherits(pseudotime, "CellDataSet"))) {
    object <- pseudotime
    pseudotime <- NULL
  }

  if (!is.null(pseudotime)) {
    if (length(pseudotime) == 1L && is.character(pseudotime) && pseudotime %in% colnames(meta)) {
      return(as.numeric(meta[[pseudotime]]))
    }
    if (length(pseudotime) == nrow(meta)) {
      return(as.numeric(pseudotime))
    }
    stop("`pseudotime` must be a metadata column name or numeric vector of length ncol(score).", call. = FALSE)
  }

  if (!is.null(object)) {
    if (inherits(object, "SlingshotDataSet") && requireNamespace("slingshot", quietly = TRUE)) {
      pt <- tryCatch(slingshot::slingPseudotime(object), error = function(e) NULL)
      if (!is.null(pt)) return(as.numeric(rowMeans(pt, na.rm = TRUE)))
    }
    if (inherits(object, "CellDataSet") && requireNamespace("monocle", quietly = TRUE)) {
      pd <- tryCatch({
        if (requireNamespace("Biobase", quietly = TRUE)) {
          Biobase::pData(object)
        } else {
          as.data.frame(object@phenoData@data)
        }
      }, error = function(e) NULL)
      if (!is.null(pd)) {
        for (nm in c("Pseudotime", "pseudotime")) {
          if (nm %in% colnames(pd)) return(as.numeric(pd[[nm]]))
        }
      }
    }
    if (inherits(object, "cell_data_set") && requireNamespace("monocle3", quietly = TRUE)) {
      pt <- tryCatch(monocle3::pseudotime(object), error = function(e) NULL)
      if (!is.null(pt)) return(as.numeric(pt))
    }
  }

  if ("pseudotime" %in% colnames(meta)) {
    return(as.numeric(meta$pseudotime))
  }

  stop("Unable to resolve pseudotime. Provide `pseudotime` explicitly (vector/column or monocle/monocle3/slingshot object).", call. = FALSE)
}

#' Extract lineage labels
#'
#' @param score `gleam_score` object.
#' @param lineage Metadata column name, vector, or NULL.
#' @param object Optional trajectory object.
#'
#' @return Character vector.
#' @export
extract_lineage <- function(score, lineage = NULL, object = NULL) {
  check_score_object(score)
  meta <- score$meta

  if (!is.null(lineage) && is.null(object) && (inherits(lineage, "SlingshotDataSet") || inherits(lineage, "cell_data_set") || inherits(lineage, "CellDataSet"))) {
    object <- lineage
    lineage <- NULL
  }

  if (!is.null(lineage)) {
    if (length(lineage) == 1L && is.character(lineage) && lineage %in% colnames(meta)) {
      return(as.character(meta[[lineage]]))
    }
    if (length(lineage) == nrow(meta)) {
      return(as.character(lineage))
    }
    stop("`lineage` must be a metadata column name or vector of length ncol(score).", call. = FALSE)
  }

  if (!is.null(object) && inherits(object, "SlingshotDataSet") && requireNamespace("slingshot", quietly = TRUE)) {
    cl <- tryCatch(slingshot::slingClusterLabels(object), error = function(e) NULL)
    if (!is.null(cl)) {
      return(apply(cl, 1, function(v) names(which.max(v))))
    }
  }
  if (!is.null(object) && inherits(object, "CellDataSet") && requireNamespace("monocle", quietly = TRUE)) {
    pd <- tryCatch({
      if (requireNamespace("Biobase", quietly = TRUE)) {
        Biobase::pData(object)
      } else {
        as.data.frame(object@phenoData@data)
      }
    }, error = function(e) NULL)
    if (!is.null(pd)) {
      for (nm in c("State", "state", "lineage", "branch")) {
        if (nm %in% colnames(pd)) return(as.character(pd[[nm]]))
      }
    }
  }

  if ("lineage" %in% colnames(meta)) {
    return(as.character(meta$lineage))
  }

  rep("lineage1", nrow(meta))
}

#' Convert trajectory inputs to a unified data.frame
#'
#' @param score `gleam_score` object.
#' @param pseudotime Pseudotime source.
#' @param lineage Lineage source.
#' @param embeddings Optional embedding matrix.
#' @param backend Trajectory backend. `auto` detects from provided inputs.
#'
#' @return Data.frame with `cell_id`, `pseudotime`, `lineage`, and embedding columns.
#' @export
as_trajectory_data <- function(
  score,
  pseudotime = NULL,
  lineage = NULL,
  embeddings = NULL,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
) {
  check_score_object(score)
  tr <- as_trajectory_result(score = score, pseudotime = pseudotime, lineage = lineage, backend = backend)

  df <- data.frame(
    cell_id = colnames(score$score),
    pseudotime = tr$pseudotime,
    lineage = tr$lineage,
    stringsAsFactors = FALSE,
    row.names = colnames(score$score)
  )

  if (!is.null(embeddings)) {
    emb <- as.data.frame(embeddings)
    if (nrow(emb) == nrow(df)) {
      emb <- emb[match(df$cell_id, rownames(emb)), , drop = FALSE]
      df <- cbind(df, emb)
    }
  }

  df
}

#' Join score matrix and trajectory metadata
#'
#' @param score `gleam_score` object.
#' @param pseudotime Pseudotime source.
#' @param lineage Lineage source.
#' @param backend Trajectory backend. `auto` detects from provided inputs.
#'
#' @return Long-format data.frame.
#' @export
join_score_trajectory <- function(score, pseudotime = NULL, lineage = NULL, backend = c("auto", "internal", "monocle3", "monocle", "slingshot")) {
  check_score_object(score)
  tr <- as_trajectory_data(score, pseudotime = pseudotime, lineage = lineage, backend = backend)
  long <- pivot_scores_long(score)
  merge(long, tr, by = "cell_id", all.x = TRUE)
}

#' Map scores to trajectory table
#'
#' @param score `gleam_score` object.
#' @param pseudotime Pseudotime source.
#' @param lineage Lineage source.
#' @param backend Trajectory backend. `auto` detects from provided inputs.
#'
#' @return Long-format data.frame.
#' @export
map_scores_to_trajectory <- function(score, pseudotime = NULL, lineage = NULL, backend = c("auto", "internal", "monocle3", "monocle", "slingshot")) {
  join_score_trajectory(score, pseudotime = pseudotime, lineage = lineage, backend = backend)
}

#' Test pathway association with trajectory
#'
#' @param score `gleam_score` object.
#' @param pathway Optional pathway name. `NULL` tests all pathways.
#' @param pseudotime Pseudotime source.
#' @param lineage Lineage source.
#' @param method One of `spearman`, `lm`, `tradeSeq`.
#' @param adjust_method P-value adjustment method.
#' @param verbose Print messages.
#' @param backend Trajectory backend. `auto` detects from provided inputs.
#'
#' @return `gleam_test` object.
#' @keywords internal
test_pathway_trajectory <- function(
  score,
  pathway = NULL,
  pseudotime = NULL,
  lineage = NULL,
  method = c("spearman", "lm", "tradeSeq"),
  adjust_method = "BH",
  verbose = TRUE,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
) {
  check_score_object(score)
  method <- match.arg(method)
  backend <- match.arg(backend)
  tr <- as_trajectory_result(score = score, pseudotime = pseudotime, lineage = lineage, backend = backend)
  pt <- tr$pseudotime
  ln <- tr$lineage
  mat <- score$score

  if (!is.null(pathway)) {
    if (!pathway %in% rownames(mat)) stop("`pathway` not found in score matrix.", call. = FALSE)
    mat <- mat[pathway, , drop = FALSE]
  }

  if (method == "tradeSeq") {
    check_method_dependency("tradeSeq", "tradeSeq")
    warning("tradeSeq interface is provided as optional placeholder in v0.2. Falling back to spearman.", call. = FALSE)
    method <- "spearman"
  }

  rows <- lapply(seq_len(nrow(mat)), function(i) {
    y <- as.numeric(mat[i, ])
    keep <- !is.na(y) & !is.na(pt)
    yk <- y[keep]
    pk <- pt[keep]

    if (length(yk) < 5) {
      return(data.frame(
        pathway = rownames(mat)[i],
        comparison_type = "trajectory",
        group1 = "pseudotime",
        group2 = NA_character_,
        celltype = NA_character_,
        level = "trajectory",
        effect_size = NA_real_,
        median_group1 = NA_real_,
        median_group2 = NA_real_,
        diff_median = NA_real_,
        p_value = NA_real_,
        p_adj = NA_real_,
        n_group1 = length(yk),
        n_group2 = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    if (method == "spearman") {
      ct <- suppressWarnings(stats::cor.test(yk, pk, method = "spearman"))
      eff <- unname(ct$estimate)
      pv <- ct$p.value
    } else {
      fit <- stats::lm(yk ~ pk)
      sm <- summary(fit)$coefficients
      eff <- sm[2, 1]
      pv <- sm[2, 4]
    }

    data.frame(
      pathway = rownames(mat)[i],
      comparison_type = "trajectory",
      group1 = "pseudotime",
      group2 = NA_character_,
      celltype = names(sort(table(as.character(ln[keep])), decreasing = TRUE))[1],
      level = "trajectory",
      effect_size = eff,
      median_group1 = stats::median(yk, na.rm = TRUE),
      median_group2 = NA_real_,
      diff_median = eff,
      p_value = pv,
      p_adj = NA_real_,
      n_group1 = length(yk),
      n_group2 = NA_integer_,
      stringsAsFactors = FALSE
    )
  })

  tbl <- do.call(rbind, rows)
  tbl$p_adj <- p_adjust_safe(tbl$p_value, adjust_method)

  if (verbose) {
    message(sprintf("[GLEAM] trajectory test completed for %d pathways.", nrow(tbl)))
  }

  new_scpathway_test(
    table = tbl,
    level = "trajectory",
    method = method,
    comparison = list(type = "trajectory", backend = tr$backend, trajectory = tr),
    params = list(adjust_method = adjust_method, backend = tr$backend, trajectory_version = tr$version_info)
  )
}
