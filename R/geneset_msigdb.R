#' Get MSigDB genesets via msigdbr
#'
#' @param species Species label accepted by `msigdbr`.
#' @param collection MSigDB collection, e.g. `H`, `C2`, `C5`.
#' @param subcollection Optional subcollection, e.g. `CP:KEGG`.
#'
#' @return Named list genesets.
#' @keywords internal
get_geneset_msigdb <- function(species = "Homo sapiens", collection = "H", subcollection = NULL) {
  require_optional_package("msigdbr", feature = "MSigDB geneset loading")

  df <- msigdbr::msigdbr(species = species, category = collection, subcategory = subcollection)
  if (nrow(df) == 0) {
    stop("No MSigDB records found for selected collection/subcollection/species.", call. = FALSE)
  }

  gs <- split(as.character(df$gene_symbol), as.character(df$gs_name))
  gs <- lapply(gs, unique)
  set_geneset_metadata(
    gs,
    source = "msigdb",
    collection = collection,
    subcollection = subcollection %||% NA_character_,
    species = species,
    version = if ("db_version" %in% colnames(df)) as.character(df$db_version[[1]]) else NA_character_
  )
}
