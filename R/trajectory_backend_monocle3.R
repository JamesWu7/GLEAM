#' @keywords internal
.monocle3_version_info <- function() {
  if (!.has_pkg("monocle3")) return(NA_character_)
  as.character(utils::packageVersion("monocle3"))
}

#' @keywords internal
.as_monocle3_input <- function(cds) {
  if (!inherits(cds, "cell_data_set")) {
    stop("`cds` must be a monocle3 `cell_data_set` object.", call. = FALSE)
  }
  cds
}

#' @keywords internal
.extract_monocle3_pseudotime <- function(cds) {
  .assert_monocle3(feature = "monocle3 trajectory backend")
  pt <- tryCatch(monocle3::pseudotime(cds), error = function(e) NULL)
  if (is.null(pt)) {
    stop("Failed to extract pseudotime from monocle3 object. Ensure pseudotime has been computed in the provided cell_data_set.", call. = FALSE)
  }
  as.numeric(pt)
}

#' @keywords internal
.extract_monocle3_lineage <- function(cds, fallback = NULL) {
  .assert_monocle3(feature = "monocle3 trajectory backend")
  cd <- tryCatch(as.data.frame(cds@colData), error = function(e) NULL)
  if (!is.null(cd)) {
    for (nm in c("lineage", "partition", "cluster", "cell_type")) {
      if (nm %in% colnames(cd)) {
        return(as.character(cd[[nm]]))
      }
    }
  }
  if (!is.null(fallback) && is.character(fallback) && length(fallback) == 1L && !is.null(cd) && fallback %in% colnames(cd)) {
    return(as.character(cd[[fallback]]))
  }
  rep("lineage1", ncol(cds))
}

#' @keywords internal
.run_trajectory_monocle3 <- function(score, cds, lineage = NULL) {
  cds <- .as_monocle3_input(cds)
  pt <- .extract_monocle3_pseudotime(cds)
  ln <- .extract_monocle3_lineage(cds, fallback = lineage)

  if (length(pt) != ncol(score$score)) {
    stop(
      "Monocle3 pseudotime length does not match score columns. Ensure score object and monocle3 object represent the same cells in the same order.",
      call. = FALSE
    )
  }

  new_gleam_trajectory_result(
    backend = "monocle3",
    pseudotime = pt,
    lineage = ln,
    meta = score$meta,
    params = list(backend = "monocle3"),
    messages = "Monocle3 backend results normalized to GLEAM trajectory object.",
    backend_object = cds,
    version_info = list(monocle3 = .monocle3_version_info())
  )
}
