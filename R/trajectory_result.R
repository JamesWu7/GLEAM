#' @keywords internal
new_gleam_trajectory_result <- function(
  backend,
  pseudotime,
  lineage,
  branch = NULL,
  meta,
  params = list(),
  messages = character(0),
  backend_object = NULL,
  version_info = list()
) {
  if (length(pseudotime) != nrow(meta)) {
    stop("`pseudotime` length must equal nrow(meta).", call. = FALSE)
  }
  if (length(lineage) != nrow(meta)) {
    stop("`lineage` length must equal nrow(meta).", call. = FALSE)
  }

  structure(
    list(
      backend = as.character(backend),
      pseudotime = as.numeric(pseudotime),
      lineage = as.character(lineage),
      branch = if (is.null(branch)) NULL else as.character(branch),
      meta = meta,
      params = params,
      messages = as.character(messages),
      backend_object = backend_object,
      version_info = version_info
    ),
    class = "gleam_trajectory_result"
  )
}

#' @keywords internal
as_trajectory_result <- function(score, pseudotime = NULL, lineage = NULL, backend = c("auto", "internal", "monocle3", "monocle", "slingshot")) {
  check_score_object(score)
  backend <- match.arg(backend)
  meta <- score$meta

  detect_backend <- function() {
    obj <- if (inherits(pseudotime, c("cell_data_set", "CellDataSet", "SlingshotDataSet"))) pseudotime else if (inherits(lineage, c("cell_data_set", "CellDataSet", "SlingshotDataSet"))) lineage else NULL
    if (is.null(obj)) return("internal")
    if (inherits(obj, "cell_data_set")) return("monocle3")
    if (inherits(obj, "CellDataSet")) return("monocle")
    if (inherits(obj, "SlingshotDataSet")) return("slingshot")
    "internal"
  }

  use_backend <- if (backend == "auto") detect_backend() else backend

  if (use_backend == "monocle3") {
    obj <- if (inherits(pseudotime, "cell_data_set")) pseudotime else if (inherits(lineage, "cell_data_set")) lineage else NULL
    if (is.null(obj)) {
      stop("Monocle3 backend selected but no monocle3 `cell_data_set` object was provided in `pseudotime` or `lineage`.", call. = FALSE)
    }
    return(.run_trajectory_monocle3(score = score, cds = obj, lineage = lineage))
  }

  pt <- extract_pseudotime(score, pseudotime = pseudotime)
  ln <- extract_lineage(score, lineage = lineage)

  if (use_backend == "monocle" && !.has_pkg("monocle")) {
    stop("Backend 'monocle' requested but optional package 'monocle' is not installed. Install with BiocManager::install('monocle').", call. = FALSE)
  }
  if (use_backend == "slingshot" && !.has_pkg("slingshot")) {
    stop("Backend 'slingshot' requested but optional package 'slingshot' is not installed. Install with BiocManager::install('slingshot').", call. = FALSE)
  }

  new_gleam_trajectory_result(
    backend = use_backend,
    pseudotime = pt,
    lineage = ln,
    meta = meta,
    params = list(backend = use_backend),
    version_info = list(
      monocle = if (.has_pkg("monocle")) as.character(utils::packageVersion("monocle")) else NA_character_,
      slingshot = if (.has_pkg("slingshot")) as.character(utils::packageVersion("slingshot")) else NA_character_
    )
  )
}
