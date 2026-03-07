#' Get geneset collection
#'
#' @param geneset Geneset input. Supported: built-in name, named list,
#'   GMT file path, or data.frame with `pathway` and `gene` columns.
#'
#' @return Named list of pathways to genes.
#' @export
get_geneset <- function(geneset) {
  if (length(geneset) == 1L && is.character(geneset) && !file.exists(geneset)) {
    if (geneset %in% c("hallmark", "immune_small")) {
      data(list = geneset, package = "scPathway", envir = environment())
      return(get(geneset, envir = environment()))
    }
  }
  as_geneset(geneset)
}

#' Read GMT file
#'
#' @param path Path to GMT file.
#'
#' @return Named list of pathways.
#' @keywords internal
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
#' @keywords internal
as_geneset <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(names(x)) || any(names(x) == "")) {
      stop("Geneset list must be a named list.", call. = FALSE)
    }
    return(lapply(x, unique))
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
#' @keywords internal
check_geneset <- function(gs, min_genes = 5, max_genes = 500) {
  if (!is.list(gs) || is.null(names(gs))) {
    stop("Geneset must be a named list.", call. = FALSE)
  }

  gs <- lapply(gs, function(g) unique(as.character(g[!is.na(g)])))
  sizes <- vapply(gs, length, integer(1))
  keep <- sizes >= min_genes & sizes <= max_genes
  if (!any(keep)) {
    stop("No pathways remained after min_genes/max_genes filtering.", call. = FALSE)
  }
  if (any(!keep)) {
    warnf("Filtered %d pathways outside gene count limits.", sum(!keep))
  }
  gs[keep]
}

#' Match genesets to expression genes
#'
#' @param gs Named list genesets.
#' @param expr_genes Character vector of expression rownames.
#' @param verbose Print overlap summary.
#'
#' @return List with matched genesets and overlap info.
#' @keywords internal
match_geneset <- function(gs, expr_genes, verbose = TRUE) {
  matched <- lapply(gs, function(g) intersect(g, expr_genes))
  n_before <- vapply(gs, length, integer(1))
  n_after <- vapply(matched, length, integer(1))
  info <- data.frame(
    pathway = names(gs),
    n_genes_input = n_before,
    n_genes_matched = n_after,
    frac_matched = n_after / pmax(1, n_before),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    message(sprintf("[scPathway] matched pathways: %d", length(gs)))
    message(sprintf("[scPathway] median matched genes: %.1f", stats::median(n_after)))
  }

  low <- n_after < 3L
  if (any(low)) {
    warnf("%d pathways have fewer than 3 matched genes.", sum(low))
  }

  list(geneset = matched, info = info)
}
