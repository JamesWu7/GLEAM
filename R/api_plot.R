#' @keywords internal
.theme_scpathway <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}
