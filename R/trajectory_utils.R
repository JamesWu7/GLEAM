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
    if (inherits(object, "cell_data_set") && requireNamespace("monocle3", quietly = TRUE)) {
      pt <- tryCatch(monocle3::pseudotime(object), error = function(e) NULL)
      if (!is.null(pt)) return(as.numeric(pt))
    }
  }

  if ("pseudotime" %in% colnames(meta)) {
    return(as.numeric(meta$pseudotime))
  }

  stop("Unable to resolve pseudotime. Provide `pseudotime` explicitly.", call. = FALSE)
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
#'
#' @return Data.frame with `cell_id`, `pseudotime`, `lineage`, and embedding columns.
#' @export
as_trajectory_data <- function(score, pseudotime = NULL, lineage = NULL, embeddings = NULL) {
  check_score_object(score)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  df <- data.frame(
    cell_id = colnames(score$score),
    pseudotime = pt,
    lineage = ln,
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
#'
#' @return Long-format data.frame.
#' @export
join_score_trajectory <- function(score, pseudotime = NULL, lineage = NULL) {
  check_score_object(score)
  tr <- as_trajectory_data(score, pseudotime = pseudotime, lineage = lineage)
  long <- pivot_scores_long(score)
  merge(long, tr, by = "cell_id", all.x = TRUE)
}

#' Map scores to trajectory table
#'
#' @param score `gleam_score` object.
#' @param pseudotime Pseudotime source.
#' @param lineage Lineage source.
#'
#' @return Long-format data.frame.
#' @export
map_scores_to_trajectory <- function(score, pseudotime = NULL, lineage = NULL) {
  join_score_trajectory(score, pseudotime = pseudotime, lineage = lineage)
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
#'
#' @return `gleam_test` object.
#' @export
test_pathway_trajectory <- function(
  score,
  pathway = NULL,
  pseudotime = NULL,
  lineage = NULL,
  method = c("spearman", "lm", "tradeSeq"),
  adjust_method = "BH",
  verbose = TRUE
) {
  check_score_object(score)
  method <- match.arg(method)
  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)
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
      celltype = as.character(stats::median(as.numeric(as.factor(ln)), na.rm = TRUE)),
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
    comparison = list(type = "trajectory"),
    params = list(adjust_method = adjust_method)
  )
}
