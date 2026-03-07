#' Build trajectory table from score object
#'
#' @inheritParams as_trajectory_data
#'
#' @return Data.frame.
#' @export
trajectory_table <- function(score, pseudotime = NULL, lineage = NULL, embeddings = NULL) {
  as_trajectory_data(score = score, pseudotime = pseudotime, lineage = lineage, embeddings = embeddings)
}
