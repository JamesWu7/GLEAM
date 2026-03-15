#' Canonical signature scoring interface
#'
#' @inheritParams score_pathway
#'
#' @section Method Guide:
#' \tabular{lll}{
#' Method \tab Best for \tab Notes \cr
#' `rank` \tab Fast robust baseline \tab Rank aggregation concept aligned with singscore \cr
#' `UCell` \tab Large single-cell data \tab Mann-Whitney U-based robust scoring \cr
#' `AUCell` \tab Rank-enrichment scoring \tab SCENIC/AUCell implementation \cr
#' `AddModuleScore` \tab Seurat-native workflows \tab Signature minus control-bin expression \cr
#' `ssGSEA` \tab Bulk-like enrichment workflows \tab GSVA package ssGSEA implementation \cr
#' `GSVA` \tab Bulk expression matrices \tab GSVA kernel-based pathway variation \cr
#' `mean` \tab Simple signatures \tab Mean expression aggregation \cr
#' `zscore` \tab Cross-signature comparability \tab Per-gene z-score then aggregate \cr
#' `robust_mean` \tab Outlier-prone data \tab Median/MAD-scaled robust aggregation \cr
#' `ensemble` \tab Consensus scoring \tab Harmonize each method (`zscore` or `rank`) before averaging \cr
#' }
#'
#' @section Method Parameters (`method_params`):
#' Supported keys include:
#' `auc_max_rank` (AUCell/native auc), `alpha` or `ssgsea_alpha` (ssGSEA),
#' `nbin`/`ctrl`/`seed` (AddModuleScore), `kcdf` (GSVA), `ucell_max_rank` (UCell).
#' For ensemble, use `ensemble_methods`, `ensemble_standardize`, and optional
#' `ensemble_weights`.
#' @return An object of class `gleam_score`.
#' @export
score_signature <- function(
  object = NULL,
  expr = NULL,
  meta = NULL,
  geneset,
  geneset_source = c("auto", "builtin", "list", "gmt", "data.frame", "msigdb", "go", "kegg", "reactome"),
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  ontology = c("BP", "MF", "CC"),
  seurat = TRUE,
  assay = NULL,
  layer = NULL,
  slot = NULL,
  method = c(
    "rank", "mean", "zscore", "scaled_mean", "robust_mean", "ensemble",
    "AddModuleScore", "UCell", "AUCell", "ssGSEA", "GSVA", "singscore",
    "auc", "zmean", "robust", "singscore_like", "ssgsea_like", "addmodulescore", "ucell", "aucell", "gsva"
  ),
  min_genes = 5,
  max_genes = 500,
  auc_max_rank = 0.05,
  ensemble_methods = c("rank", "zscore", "mean"),
  ensemble_combine = c("mean", "median"),
  ensemble_standardize = c("zscore", "rank"),
  ensemble_weights = NULL,
  method_params = list(),
  verbose = TRUE
) {
  score_pathway(
    object = object,
    expr = expr,
    meta = meta,
    geneset = geneset,
    geneset_source = geneset_source,
    species = species,
    collection = collection,
    subcollection = subcollection,
    ontology = ontology,
    seurat = seurat,
    assay = assay,
    layer = layer,
    slot = slot,
    method = method,
    min_genes = min_genes,
    max_genes = max_genes,
    auc_max_rank = auc_max_rank,
    ensemble_methods = ensemble_methods,
    ensemble_combine = ensemble_combine,
    ensemble_standardize = ensemble_standardize,
    ensemble_weights = ensemble_weights,
    method_params = method_params,
    verbose = verbose
  )
}

#' Canonical signature aggregation interface
#'
#' @inheritParams aggregate_pathway
#' @return Aggregated table (long) or signature-by-group matrix (wide).
#' @export
aggregate_signature <- function(score, by, fun = c("mean", "median", "fraction", "sum"), threshold = 0, long = TRUE) {
  aggregate_pathway(score = score, by = by, fun = fun, threshold = threshold, long = long)
}

#' Canonical signature testing interface
#'
#' @inheritParams test_pathway
#' @return An object of class `gleam_test`.
#' @export
test_signature <- function(
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
  backend <- match.arg(backend)
  test_pathway(
    score = score,
    group = group,
    sample = sample,
    celltype = celltype,
    target_celltype = target_celltype,
    region = region,
    pseudotime = pseudotime,
    lineage = lineage,
    level = level,
    method = method,
    ref_group = ref_group,
    paired = paired,
    covariates = covariates,
    adjust_method = adjust_method,
    min_cells = min_cells,
    min_samples = min_samples,
    aggregation = aggregation,
    threshold = threshold,
    backend = backend,
    verbose = verbose
  )
}

#' Canonical spatial signature testing interface
#'
#' @inheritParams test_pathway_spatial
#' @return `gleam_test` object.
#' @export
test_signature_spatial <- function(
  score,
  region,
  group = NULL,
  sample = NULL,
  method = c("wilcox", "t"),
  adjust_method = "BH",
  level = c("region", "sample_region")
) {
  test_pathway_spatial(
    score = score,
    region = region,
    group = group,
    sample = sample,
    method = method,
    adjust_method = adjust_method,
    level = level
  )
}

#' Canonical trajectory signature testing interface
#'
#' @inheritParams test_pathway_trajectory
#' @param backend Trajectory backend to use. `auto` detects from provided inputs.
#' @return `gleam_test` object.
#' @export
test_signature_trajectory <- function(
  score,
  pathway = NULL,
  pseudotime = NULL,
  lineage = NULL,
  method = c("spearman", "lm", "tradeSeq"),
  adjust_method = "BH",
  verbose = TRUE,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
) {
  backend <- match.arg(backend)
  test_pathway_trajectory(
    score = score,
    pathway = pathway,
    pseudotime = pseudotime,
    lineage = lineage,
    method = method,
    adjust_method = adjust_method,
    verbose = verbose,
    backend = backend
  )
}
