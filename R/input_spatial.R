#' Build spatial table from score object
#'
#' @param score `gleam_score` object.
#' @param coords Coordinate data.frame with x/y.
#' @param meta Optional metadata table.
#'
#' @return Data.frame.
#' @keywords internal
spatial_table <- function(score, coords, meta = NULL) {
  join_score_spatial(score = score, coords = coords, meta = meta)
}
