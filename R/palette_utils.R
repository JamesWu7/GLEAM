#' List available GLEAM palettes
#'
#' @return Character vector of palette names.
#' @export
list_palettes <- function() {
  c(
    "gleam_discrete",
    "gleam_continuous",
    "viridis",
    "brewer_set2",
    "brewer_dark2"
  )
}

#' Get palette colors
#'
#' @param name Palette name.
#' @param n Number of colors.
#' @param continuous Whether to return continuous palette function output.
#' @param reverse Reverse color order.
#'
#' @return Character color vector.
#' @export
get_palette <- function(name = "gleam_discrete", n = 8, continuous = FALSE, reverse = FALSE) {
  cols <- switch(
    name,
    gleam_discrete = {
      base <- c("#1f3b73", "#2a9d8f", "#e9c46a", "#e76f51", "#8ab17d", "#264653", "#f4a261", "#6d597a")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    gleam_continuous = grDevices::colorRampPalette(c("#1f3b73", "#f4f1de", "#e63946"))(max(2, n)),
    viridis = {
      check_plot_dependency("palette viridis", "viridisLite")
      viridisLite::viridis(max(2, n))
    },
    brewer_set2 = {
      check_plot_dependency("palette brewer_set2", "RColorBrewer")
      RColorBrewer::brewer.pal(max(3, min(8, n)), "Set2")
    },
    brewer_dark2 = {
      check_plot_dependency("palette brewer_dark2", "RColorBrewer")
      RColorBrewer::brewer.pal(max(3, min(8, n)), "Dark2")
    },
    stop(sprintf("Unknown palette '%s'. Use list_palettes().", name), call. = FALSE)
  )

  if (continuous && length(cols) < n) {
    cols <- grDevices::colorRampPalette(cols)(n)
  }
  if (reverse) cols <- rev(cols)
  cols
}

#' GLEAM color scale helper
#'
#' @param palette Palette name or custom color vector.
#' @param continuous Whether to return continuous scale.
#'
#' @return A ggplot2 scale object.
#' @export
scale_gleam_color <- function(palette = "gleam_discrete", continuous = FALSE) {
  cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 9, continuous = continuous) else palette
  if (continuous) {
    ggplot2::scale_color_gradientn(colors = cols)
  } else {
    ggplot2::scale_color_manual(values = cols)
  }
}

#' GLEAM fill scale helper
#'
#' @inheritParams scale_gleam_color
#' @return A ggplot2 scale object.
#' @export
scale_gleam_fill <- function(palette = "gleam_discrete", continuous = FALSE) {
  cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 9, continuous = continuous) else palette
  if (continuous) {
    ggplot2::scale_fill_gradientn(colors = cols)
  } else {
    ggplot2::scale_fill_manual(values = cols)
  }
}
