#' Build trajectory table from score object
#'
#' @inheritParams as_trajectory_data
#'
#' @return Data.frame.
#' @keywords internal
trajectory_table <- function(
  score,
  pseudotime = NULL,
  lineage = NULL,
  embeddings = NULL,
  backend = c("auto", "internal", "monocle3", "monocle", "slingshot")
) {
  as_trajectory_data(
    score = score,
    pseudotime = pseudotime,
    lineage = lineage,
    embeddings = embeddings,
    backend = backend
  )
}
