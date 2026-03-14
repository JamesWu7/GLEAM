#' GLEAM plotting theme
#'
#' @param base_size Base text size used for the overall theme.
#' @param title_size Plot title size.
#' @param subtitle_size Plot subtitle size.
#' @param axis_title_size Axis title size.
#' @param axis_text_size Axis text size.
#' @param legend_title_size Legend title size.
#' @param legend_text_size Legend text size.
#' @param strip_text_size Facet strip text size.
#' @param font_family Font family.
#' @param font_face Font face.
#' @param title_color Title text color.
#' @param subtitle_color Subtitle text color.
#' @param axis_text_color Axis text color.
#' @param axis_title_color Axis title color.
#' @param legend_text_color Legend text color.
#' @param legend_title_color Legend title color.
#' @param strip_text_color Facet strip text color.
#' @param text_color Global text color fallback.
#'
#' @return A ggplot2 theme object.
#' @export
gleam_theme <- function(
  base_size = 14,
  title_size = base_size + 4,
  subtitle_size = base_size + 1,
  axis_title_size = base_size + 1,
  axis_text_size = max(10, base_size),
  legend_title_size = base_size + 1,
  legend_text_size = max(10, base_size - 1),
  strip_text_size = base_size + 1,
  font_family = "",
  font_face = "plain",
  title_color = "#1F2937",
  subtitle_color = "#374151",
  axis_text_color = "#111827",
  axis_title_color = "#1F2937",
  legend_text_color = "#111827",
  legend_title_color = "#1F2937",
  strip_text_color = "#1F2937",
  text_color = "#111827"
) {
  ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(color = text_color, face = font_face, family = font_family),
      plot.title = ggplot2::element_text(size = title_size, color = title_color, face = "bold", family = font_family),
      plot.subtitle = ggplot2::element_text(size = subtitle_size, color = subtitle_color, family = font_family),
      axis.title = ggplot2::element_text(size = axis_title_size, color = axis_title_color, face = "bold", family = font_family),
      axis.text = ggplot2::element_text(size = axis_text_size, color = axis_text_color, family = font_family),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.title = ggplot2::element_text(size = legend_title_size, color = legend_title_color, face = "bold", family = font_family),
      legend.text = ggplot2::element_text(size = legend_text_size, color = legend_text_color, family = font_family),
      strip.text = ggplot2::element_text(size = strip_text_size, color = strip_text_color, face = "bold", family = font_family),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Apply GLEAM theme to a ggplot
#'
#' @param p A ggplot object.
#' @param ... Parameters passed to [gleam_theme()].
#'
#' @return A ggplot object.
#' @export
apply_gleam_theme <- function(p, ...) {
  p + gleam_theme(...)
}

#' @keywords internal
resolve_text_params <- function(x = list()) {
  if (is.null(x)) return(list())
  if (!is.list(x)) stop("`theme_params` must be a list.", call. = FALSE)
  x
}

#' @keywords internal
resolve_palette_params <- function(palette = NULL, fill_palette = NULL, color_palette = NULL) {
  list(
    palette = palette %||% "gleam_discrete",
    fill_palette = fill_palette %||% palette %||% "gleam_discrete",
    color_palette = color_palette %||% palette %||% "gleam_continuous"
  )
}

#' @keywords internal
.theme_gleam <- function(base_size = 13) {
  gleam_theme(base_size = base_size)
}
