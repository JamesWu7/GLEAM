#' @keywords internal
get_geneset_go <- function(species = "Homo sapiens", ontology = c("BP", "MF", "CC")) {
  ontology <- match.arg(ontology)
  get_geneset_msigdb(
    species = species,
    collection = "C5",
    subcollection = paste0("GO:", ontology)
  )
}

#' @keywords internal
get_geneset_kegg <- function(species = "Homo sapiens") {
  get_geneset_msigdb(
    species = species,
    collection = "C2",
    subcollection = "CP:KEGG"
  )
}

#' @keywords internal
get_geneset_reactome <- function(species = "Homo sapiens") {
  get_geneset_msigdb(
    species = species,
    collection = "C2",
    subcollection = "CP:REACTOME"
  )
}
