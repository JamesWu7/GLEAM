#' Built-in hallmark-like gene sets
#'
#' A named list of hallmark-like pathways and genes.
#'
#' @format A named list.
#' @source Internal curated collection
"hallmark"

#' Built-in immune-focused small gene sets
#'
#' A named list of immune pathways and genes.
#'
#' @format A named list.
#' @source Internal curated collection
"immune_small"

#' Tiny toy expression dataset for tests
#'
#' A list with expression matrix and metadata for fast unit tests.
#'
#' @format A list with components `expr` and `meta`.
#' @source Simulated data
"toy_expr"

#' Small PBMC example dataset (legacy compatibility)
#'
#' A list with sparse expression matrix and metadata.
#'
#' @format A list with components `expr` and `meta`.
#' @source Simulated PBMC-like data
"pbmc_small"

#' Medium PBMC expression matrix for realistic demos
#'
#' Sparse gene-by-cell matrix suitable for tutorial workflows.
#'
#' @format A `dgCMatrix`.
#' @source Simulated PBMC-like data
"pbmc_medium_matrix"

#' Medium PBMC metadata for realistic demos
#'
#' Metadata aligned to `pbmc_medium_matrix` columns.
#'
#' @format A `data.frame`.
#' @source Simulated PBMC-like data
"pbmc_medium_meta"

#' Medium spatial expression matrix
#'
#' Sparse gene-by-spot matrix for spatial workflow demos.
#'
#' @format A `dgCMatrix`.
#' @source Simulated spatial transcriptomics-like data
"spatial_medium_expr"

#' Medium spatial metadata
#'
#' Metadata for spatial spots including sample/region/domain and coordinates.
#'
#' @format A `data.frame`.
#' @source Simulated spatial transcriptomics-like data
"spatial_medium_meta"

#' Medium spatial coordinates
#'
#' Coordinates aligned to spatial spots.
#'
#' @format A `data.frame` with columns `x` and `y`.
#' @source Simulated spatial transcriptomics-like data
"spatial_medium_coords"
