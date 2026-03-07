#' List geneset sources
#'
#' @return Data.frame of supported geneset sources.
#' @export
list_geneset_sources <- function() {
  data.frame(
    source = c("builtin", "list", "gmt", "data.frame", "msigdb", "go", "kegg", "reactome"),
    requires_package = c(NA, NA, NA, NA, "msigdbr", "msigdbr", "msigdbr", "msigdbr"),
    supported_builtin_species = c("human,mouse", "any(custom)", "any(custom)", "any(custom)", "human,mouse", "human,mouse", "human,mouse", "human,mouse"),
    built_in_available = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    internet_required = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

#' Search pathways in a geneset collection
#'
#' @param geneset Geneset input accepted by [get_geneset()].
#' @param pattern Search pattern.
#' @param ignore.case Whether search is case-insensitive.
#' @param ... Additional arguments passed to [get_geneset()].
#'
#' @return Character vector of matched pathway names.
#' @export
search_geneset <- function(geneset = "hallmark", pattern, ignore.case = TRUE, ...) {
  gs <- get_geneset(geneset = geneset, ...)
  names(gs)[grepl(pattern, names(gs), ignore.case = ignore.case)]
}
