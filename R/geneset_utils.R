#' Get geneset collection
#'
#' @param geneset Geneset input.
#' @param source Geneset source. One of `auto`, `builtin`, `list`, `gmt`,
#'   `data.frame`, `msigdb`, `go`, `kegg`, `reactome`.
#' @param species Species label used by external geneset sources.
#' @param collection MSigDB collection (for source = `msigdb`).
#' @param subcollection MSigDB subcollection (for source = `msigdb`).
#' @param ontology GO ontology (`BP`, `MF`, `CC`) for source = `go`.
#' @param version Optional source version string.
#'
#' @return Named list of pathways to genes with source metadata attached.
#' @export
get_geneset <- function(
  geneset,
  source = c("auto", "builtin", "list", "gmt", "data.frame", "msigdb", "go", "kegg", "reactome"),
  species = "Homo sapiens",
  collection = "H",
  subcollection = NULL,
  ontology = c("BP", "MF", "CC"),
  version = NA_character_
) {
  source <- match.arg(source)
  ontology <- match.arg(ontology)

  if (source == "auto") {
    source <- infer_geneset_source(geneset)
  }

  gs <- switch(
    source,
    builtin = get_builtin_geneset(geneset),
    list = as_geneset(geneset),
    gmt = read_gmt(as.character(geneset)),
    `data.frame` = as_geneset(geneset),
    msigdb = get_geneset_msigdb(species = species, collection = collection, subcollection = subcollection),
    go = get_geneset_go(species = species, ontology = ontology),
    kegg = get_geneset_kegg(species = species),
    reactome = get_geneset_reactome(species = species),
    stop(sprintf("Unsupported geneset source: %s", source), call. = FALSE)
  )

  set_geneset_metadata(
    gs,
    source = source,
    collection = collection,
    subcollection = subcollection,
    species = species,
    version = version
  )
}

#' Read GMT file
#'
#' @param path Path to GMT file.
#'
#' @return Named list of pathways.
#' @export
read_gmt <- function(path) {
  if (!file.exists(path)) {
    stop("GMT file does not exist: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  out <- vector("list", length(lines))
  nms <- character(length(lines))
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[[i]], "\t", fixed = FALSE)[[1]]
    if (length(parts) < 3L) {
      next
    }
    nms[[i]] <- parts[[1]]
    out[[i]] <- unique(parts[-c(1, 2)])
  }
  keep <- nzchar(nms)
  out <- out[keep]
  nms <- nms[keep]
  names(out) <- nms
  out
}

#' Coerce geneset input to named list
#'
#' @param x Geneset input.
#'
#' @return Named list.
#' @export
as_geneset <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("Geneset list must be a named list.", call. = FALSE)
    }
    return(lapply(x, function(v) unique(as.character(v[!is.na(v)]))))
  }

  if (length(x) == 1L && is.character(x) && file.exists(x)) {
    return(read_gmt(x))
  }

  if (is.data.frame(x)) {
    req <- c("pathway", "gene")
    miss <- setdiff(req, colnames(x))
    if (length(miss) > 0L) {
      stop("Geneset data.frame must include columns: pathway, gene.", call. = FALSE)
    }
    sp <- split(as.character(x$gene), as.character(x$pathway))
    return(lapply(sp, unique))
  }

  stop("Unsupported geneset input. Use built-in name, named list, GMT path, or data.frame.", call. = FALSE)
}

#' Validate genesets by size
#'
#' @param gs Named list genesets.
#' @param min_genes Minimum genes per pathway.
#' @param max_genes Maximum genes per pathway.
#'
#' @return Filtered named list.
#' @export
check_geneset <- function(gs, min_genes = 5, max_genes = 500) {
  if (!is.list(gs) || is.null(names(gs))) {
    stop("Geneset must be a named list.", call. = FALSE)
  }

  meta <- get_geneset_metadata(gs)
  gs <- lapply(gs, function(g) unique(as.character(g[!is.na(g)])))
  sizes <- vapply(gs, length, integer(1))
  keep <- sizes >= min_genes & sizes <= max_genes
  if (!any(keep)) {
    stop("No pathways remained after min_genes/max_genes filtering.", call. = FALSE)
  }
  if (any(!keep)) {
    warnf("Filtered %d pathways outside gene count limits.", sum(!keep))
  }
  gs <- gs[keep]
  if (!is.null(meta)) {
    meta <- meta[meta$pathway %in% names(gs), , drop = FALSE]
    gs <- set_geneset_metadata(gs, metadata = meta)
  }
  gs
}

#' Match genesets to expression genes
#'
#' @param gs Named list genesets.
#' @param expr_genes Character vector of expression rownames.
#' @param verbose Print overlap summary.
#'
#' @return List with matched genesets and overlap info.
#' @export
match_geneset <- function(gs, expr_genes, verbose = TRUE) {
  matched <- lapply(gs, function(g) intersect(g, expr_genes))
  n_before <- vapply(gs, length, integer(1))
  n_after <- vapply(matched, length, integer(1))

  base_meta <- get_geneset_metadata(gs)
  info <- data.frame(
    pathway = names(gs),
    n_genes_input = n_before,
    n_genes_matched = n_after,
    frac_matched = n_after / pmax(1, n_before),
    stringsAsFactors = FALSE
  )
  if (!is.null(base_meta)) {
    info <- merge(base_meta, info, by = "pathway", all.y = TRUE, sort = FALSE)
  }

  if (verbose) {
    message(sprintf("[GLEAM] matched pathways: %d", length(gs)))
    message(sprintf("[GLEAM] median matched genes: %.1f", stats::median(n_after)))
  }

  low <- n_after < 3L
  if (any(low)) {
    warnf("%d pathways have fewer than 3 matched genes.", sum(low))
  }

  matched <- set_geneset_metadata(matched, metadata = info)
  list(geneset = matched, info = info)
}

#' @keywords internal
infer_geneset_source <- function(geneset) {
  if (length(geneset) == 1L && is.character(geneset) && geneset %in% c("hallmark", "immune_small")) return("builtin")
  if (is.list(geneset) && !is.data.frame(geneset)) return("list")
  if (length(geneset) == 1L && is.character(geneset) && file.exists(geneset)) return("gmt")
  if (is.data.frame(geneset)) return("data.frame")
  "list"
}

#' @keywords internal
get_builtin_geneset <- function(geneset) {
  if (!(length(geneset) == 1L && is.character(geneset))) {
    stop("Built-in geneset must be provided as a character name.", call. = FALSE)
  }
  if (!geneset %in% c("hallmark", "immune_small")) {
    stop(sprintf("Unknown built-in geneset '%s'.", geneset), call. = FALSE)
  }
  data(list = geneset, package = "GLEAM", envir = environment())
  get(geneset, envir = environment())
}

#' @keywords internal
set_geneset_metadata <- function(gs, source = NA_character_, collection = NA_character_, subcollection = NA_character_, species = NA_character_, version = NA_character_, metadata = NULL) {
  if (is.null(metadata)) {
    n <- length(gs)
    rep_len_chr <- function(x) {
      if (is.null(x) || length(x) == 0L) x <- NA_character_
      rep(as.character(x[[1]]), n)
    }
    metadata <- data.frame(
      pathway = names(gs),
      source = rep_len_chr(source),
      collection = rep_len_chr(collection),
      subcollection = rep_len_chr(subcollection),
      species = rep_len_chr(species),
      version = rep_len_chr(version),
      stringsAsFactors = FALSE
    )
  }
  attr(gs, "gleam_geneset_metadata") <- metadata
  gs
}

#' @keywords internal
get_geneset_metadata <- function(gs) {
  attr(gs, "gleam_geneset_metadata")
}
