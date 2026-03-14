# Internal helpers used by vignettes/scripts to avoid redefining ad-hoc functions.

#' @keywords internal
.pick_first_col <- function(candidates, cols) {
  hit <- candidates[candidates %in% cols]
  if (length(hit) == 0L) return(NULL)
  hit[[1]]
}

#' @keywords internal
.stratified_keep <- function(meta, n_target, strata_candidates = character()) {
  if (nrow(meta) <= n_target) return(rownames(meta))
  strata <- intersect(strata_candidates, colnames(meta))
  all_ids <- rownames(meta)
  if (length(strata) == 0L) return(all_ids[seq_len(n_target)])

  key <- interaction(meta[, strata, drop = FALSE], drop = TRUE, lex.order = TRUE)
  groups <- split(all_ids, key)
  per_group <- max(1L, floor(n_target / max(1L, length(groups))))
  keep <- unlist(lapply(groups, function(ids) {
    ids <- sort(ids)
    head(ids, per_group)
  }), use.names = FALSE)
  if (length(keep) < n_target) {
    keep <- c(keep, head(setdiff(all_ids, keep), n_target - length(keep)))
  }
  unique(keep)[seq_len(min(n_target, length(unique(keep))))]
}

#' @keywords internal
.map_geneset_to_expr <- function(gs, expr_genes, min_genes = 1L) {
  expr_genes <- as.character(expr_genes)
  expr_upper <- toupper(expr_genes)
  mapped <- lapply(gs, function(g) {
    idx <- match(toupper(unique(as.character(g))), expr_upper, nomatch = 0L)
    unique(expr_genes[idx[idx > 0L]])
  })
  mapped[vapply(mapped, length, integer(1)) >= as.integer(min_genes)]
}

#' @keywords internal
.fallback_signatures <- function(expr_genes, max_genes = 30L, min_genes = 3L) {
  g <- unique(as.character(expr_genes))
  g <- g[!is.na(g) & nzchar(g)]
  k <- min(as.integer(max_genes), max(as.integer(min_genes), floor(length(g) / 4)))
  idx1 <- seq_len(min(k, length(g)))
  idx2 <- seq.int(max(1L, length(g) - k + 1L), length(g))
  list(
    Signature_A = g[idx1],
    Signature_B = g[idx2]
  )
}

#' @keywords internal
.make_discrete_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(9, "Set1")
    if (n <= length(base)) return(base[seq_len(n)])
    return(grDevices::colorRampPalette(base)(n))
  }
  get_palette("gleam_discrete", n = n)
}

#' @keywords internal
.safe_print_plot <- function(label, expr) {
  p <- tryCatch(eval.parent(substitute(expr)), error = function(e) {
    message("[GLEAM] ", label, " skipped: ", conditionMessage(e))
    NULL
  })
  if (is.null(p)) return(invisible(NULL))
  tryCatch(print(p), error = function(e) {
    message("[GLEAM] ", label, " render skipped: ", conditionMessage(e))
  })
  invisible(NULL)
}

#' @keywords internal
.resolve_seurat_expr_matrix <- function(obj) {
  assay_candidates <- unique(c(
    tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) NULL),
    "Spatial",
    "RNA"
  ))
  assay_candidates <- assay_candidates[!is.na(assay_candidates) & nzchar(assay_candidates)]

  get_if_nonempty <- function(expr) {
    m <- tryCatch(eval.parent(substitute(expr)), error = function(e) NULL)
    if (!is.null(m) && nrow(m) > 0L && ncol(m) > 0L) return(m)
    NULL
  }

  for (assay in assay_candidates) {
    m <- get_if_nonempty(SeuratObject::LayerData(object = obj, assay = assay, layer = "data"))
    if (!is.null(m)) return(m)
    m <- get_if_nonempty(SeuratObject::LayerData(object = obj, assay = assay, layer = "counts"))
    if (!is.null(m)) return(m)
    m <- get_if_nonempty(SeuratObject::GetAssayData(object = obj, assay = assay, slot = "data"))
    if (!is.null(m)) return(m)
    m <- get_if_nonempty(SeuratObject::GetAssayData(object = obj, assay = assay, slot = "counts"))
    if (!is.null(m)) return(m)
  }

  stop("Failed to extract non-empty expression matrix from Seurat spatial object.")
}

#' @keywords internal
.prep_spatial_object_matrix <- function(obj, sample_label) {
  expr <- .resolve_seurat_expr_matrix(obj)
  md <- as.data.frame(obj[[]], stringsAsFactors = FALSE)
  if (is.null(rownames(md))) rownames(md) <- colnames(expr)
  md <- md[colnames(expr), , drop = FALSE]
  if (!"sample" %in% colnames(md)) md$sample <- sample_label
  if (!"section" %in% colnames(md)) md$section <- sample_label
  if (!all(c("x", "y") %in% colnames(md))) {
    if (all(c("imagecol", "imagerow") %in% colnames(md))) {
      md$x <- md$imagecol
      md$y <- md$imagerow
    } else if (all(c("col", "row") %in% colnames(md))) {
      md$x <- md$col
      md$y <- md$row
    }
  }
  prefixed_ids <- make.unique(paste0(sample_label, "_", colnames(expr)), sep = "_dup")
  colnames(expr) <- prefixed_ids
  rownames(md) <- prefixed_ids
  list(expr = expr, meta = md)
}
