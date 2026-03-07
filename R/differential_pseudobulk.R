#' @keywords internal
test_pathway_pseudobulk <- function(
  score,
  group,
  sample,
  celltype = NULL,
  region = NULL,
  method = c("wilcox", "t", "lm", "limma", "edgeR", "DESeq2"),
  adjust_method = "BH",
  aggregation = c("mean", "median", "fraction", "sum"),
  threshold = 0,
  ref_group = NULL,
  level = "pseudobulk"
) {
  method <- match.arg(method)
  aggregation <- match.arg(aggregation)

  by_cols <- c("sample", "group")
  meta <- score$meta
  meta$.group <- as.character(group)
  meta$.sample <- as.character(sample)

  if (!is.null(celltype)) {
    meta$.celltype <- as.character(celltype)
    by_cols <- c(by_cols, "celltype")
  }
  if (!is.null(region)) {
    meta$.region <- as.character(region)
    by_cols <- c(by_cols, "region")
  }

  # Build transient score object with remapped metadata names expected by aggregate_pathway.
  tmp <- score
  tmp$meta <- data.frame(
    sample = meta$.sample,
    group = meta$.group,
    celltype = meta$.celltype %||% NA_character_,
    region = meta$.region %||% NA_character_,
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )

  agg <- aggregate_pathway(tmp, by = by_cols, fun = aggregation, threshold = threshold, long = TRUE)

  if (method %in% c("limma", "edgeR", "DESeq2")) {
    pkg <- switch(method, limma = "limma", edgeR = "edgeR", DESeq2 = "DESeq2")
    check_method_dependency(method, pkg)
    warning(sprintf("%s pseudobulk backend currently routes to lm-style score testing in v0.2.", method), call. = FALSE)
    method <- "lm"
  }

  rows <- lapply(split(agg, agg$pathway), function(dfp) {
    g <- as.character(dfp$group)
    if (length(unique(g)) < 2L) {
      return(data.frame(
        pathway = unique(dfp$pathway),
        comparison_type = "pseudobulk",
        group1 = NA_character_,
        group2 = NA_character_,
        celltype = if ("celltype" %in% colnames(dfp)) paste(unique(dfp$celltype), collapse = ";") else NA_character_,
        level = level,
        effect_size = NA_real_,
        median_group1 = NA_real_,
        median_group2 = NA_real_,
        diff_median = NA_real_,
        p_value = NA_real_,
        p_adj = NA_real_,
        n_group1 = NA_integer_,
        n_group2 = NA_integer_,
        mean_group1 = NA_real_,
        mean_group2 = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    gp <- pick_two_groups(g, ref_group = ref_group)
    keep <- g %in% gp
    x <- dfp$value[keep]
    g2 <- g[keep]

    if (method == "lm") {
      fit <- stats::lm(x ~ as.factor(g2))
      sm <- summary(fit)$coefficients
      pv <- if (nrow(sm) >= 2L) sm[2, 4] else NA_real_
    } else {
      pv <- run_group_test(x, g2, method = ifelse(method == "t", "t", "wilcox"))$p_value
    }

    eff <- compute_effects(x, g2)
    data.frame(
      pathway = unique(dfp$pathway),
      comparison_type = "pseudobulk",
      group1 = eff$group1,
      group2 = eff$group2,
      celltype = if ("celltype" %in% colnames(dfp)) paste(unique(dfp$celltype), collapse = ";") else NA_character_,
      level = level,
      effect_size = eff$effect_size,
      median_group1 = eff$median_group1,
      median_group2 = eff$median_group2,
      diff_median = eff$diff_median,
      p_value = pv,
      p_adj = NA_real_,
      n_group1 = sum(g2 == eff$group1, na.rm = TRUE),
      n_group2 = sum(g2 == eff$group2, na.rm = TRUE),
      mean_group1 = eff$mean_group1,
      mean_group2 = eff$mean_group2,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$p_adj <- p_adjust_safe(out$p_value, method = adjust_method)
  out$direction <- ifelse(is.na(out$effect_size), NA_character_, ifelse(out$effect_size >= 0, "up", "down"))
  rownames(out) <- NULL
  out
}
