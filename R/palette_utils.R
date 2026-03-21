#' List available GLEAM palettes
#'
#' @return Character vector of palette names.
#' @export
list_palettes <- function() {
  c(
    "gleam_discrete",
    "gleam_continuous",
    "okabe_ito",
    "tableau10",
    "viridis",
    "magma",
    "plasma",
    "brewer_set2",
    "brewer_dark2",
    "manual"
  )
}

.is_single_color_literal <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(FALSE)
  isTRUE(tryCatch({
    grDevices::col2rgb(x)
    TRUE
  }, error = function(e) FALSE))
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
  if (length(name) == 1L && .is_single_color_literal(name)) {
    cols <- rep(as.character(name), max(1L, n))
    if (reverse) cols <- rev(cols)
    return(cols)
  }

  if (length(name) > 1L) {
    cols <- as.character(name)
    if (continuous && length(cols) < n) cols <- grDevices::colorRampPalette(cols)(n)
    if (reverse) cols <- rev(cols)
    return(cols)
  }

  cols <- switch(
    name,
    gleam_discrete = {
      base <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    okabe_ito = {
      base <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    tableau10 = {
      base <- c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2", "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    gleam_continuous = grDevices::colorRampPalette(c("#1f3b73", "#f4f1de", "#e63946"))(max(2, n)),
    viridis = {
      check_plot_dependency("palette viridis", "viridisLite")
      viridisLite::viridis(max(2, n))
    },
    magma = {
      check_plot_dependency("palette magma", "viridisLite")
      viridisLite::magma(max(2, n))
    },
    plasma = {
      check_plot_dependency("palette plasma", "viridisLite")
      viridisLite::plasma(max(2, n))
    },
    brewer_set2 = {
      check_plot_dependency("palette brewer_set2", "RColorBrewer")
      base <- RColorBrewer::brewer.pal(8, "Set2")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    brewer_dark2 = {
      check_plot_dependency("palette brewer_dark2", "RColorBrewer")
      base <- RColorBrewer::brewer.pal(8, "Dark2")
      if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
    },
    {
      if (requireNamespace("paletteer", quietly = TRUE)) {
        pal <- tryCatch(as.character(paletteer::paletteer_d(name, n = n)), error = function(e) NULL)
        if (!is.null(pal)) {
          pal
        } else {
          stop(sprintf("Unknown palette '%s'. Use list_palettes() or provide a custom color vector.", name), call. = FALSE)
        }
      } else {
        stop(sprintf("Unknown palette '%s'. Use list_palettes() or install 'paletteer'.", name), call. = FALSE)
      }
    }
  )

  if (continuous && length(cols) < n) {
    cols <- grDevices::colorRampPalette(cols)(n)
  }
  if (reverse) cols <- rev(cols)
  cols
}

#' @keywords internal
resolve_discrete_levels <- function(values, palette = "gleam_discrete") {
  vals <- as.character(values)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0L) return(character(0))
  base_levels <- if (is.factor(values)) levels(base::droplevels(values)) else unique(vals)
  if (!is.character(palette) && !is.null(names(palette)) && any(nzchar(names(palette)))) {
    pal_levels <- names(palette)[nzchar(names(palette))]
    return(c(pal_levels[pal_levels %in% vals], setdiff(base_levels, pal_levels)))
  }
  base_levels
}

#' @keywords internal
resolve_discrete_palette_values <- function(levels, palette = "gleam_discrete", reverse = FALSE) {
  lv <- as.character(levels)
  lv <- lv[nzchar(lv)]
  if (length(lv) == 0L) return(stats::setNames(character(0), character(0)))
  if (!is.character(palette) && !is.null(names(palette)) && any(nzchar(names(palette)))) {
    vals <- as.character(palette)
    names(vals) <- names(palette)
    out <- vals[lv]
    miss <- is.na(out) | !nzchar(out)
    if (any(miss)) {
      out[miss] <- get_palette(name = "gleam_discrete", n = sum(miss), continuous = FALSE)
    }
    return(stats::setNames(as.character(out), lv))
  }
  cols <- if (is.character(palette) && length(palette) == 1L) {
    get_palette(name = palette, n = length(lv), continuous = FALSE, reverse = reverse)
  } else {
    pal <- as.character(palette)
    if (reverse) pal <- rev(pal)
    if (length(pal) >= length(lv)) pal[seq_len(length(lv))] else grDevices::colorRampPalette(pal)(length(lv))
  }
  stats::setNames(cols, lv)
}

#' GLEAM color scale helper
#'
#' @param palette Palette name or custom color vector.
#' @param continuous Whether to return continuous scale.
#'
#' @return A ggplot2 scale object.
#' @export
scale_gleam_color <- function(palette = "gleam_discrete", continuous = FALSE) {
  if (continuous) {
    cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 256, continuous = TRUE) else as.character(palette)
    ggplot2::scale_color_gradientn(colors = cols)
  } else {
    is_named_single <- is.character(palette) && length(palette) == 1L &&
      !is.null(names(palette)) && any(nzchar(names(palette)))
    use_named_palette <- is.character(palette) && length(palette) == 1L &&
      !is_named_single && !.is_single_color_literal(palette)
    if (use_named_palette) {
      ggplot2::discrete_scale(
        aesthetics = "colour",
        palette = function(n) get_palette(name = palette, n = n, continuous = FALSE)
      )
    } else {
      pal <- as.character(palette)
      if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
        ggplot2::scale_color_manual(values = palette, drop = FALSE)
      } else {
        ggplot2::discrete_scale(
          aesthetics = "colour",
          palette = function(n) {
            if (length(pal) >= n) pal[seq_len(n)] else grDevices::colorRampPalette(pal)(n)
          }
        )
      }
    }
  }
}

#' GLEAM fill scale helper
#'
#' @inheritParams scale_gleam_color
#' @return A ggplot2 scale object.
#' @export
scale_gleam_fill <- function(palette = "gleam_discrete", continuous = FALSE) {
  if (continuous) {
    cols <- if (is.character(palette) && length(palette) == 1L) get_palette(palette, n = 256, continuous = TRUE) else as.character(palette)
    ggplot2::scale_fill_gradientn(colors = cols)
  } else {
    is_named_single <- is.character(palette) && length(palette) == 1L &&
      !is.null(names(palette)) && any(nzchar(names(palette)))
    use_named_palette <- is.character(palette) && length(palette) == 1L &&
      !is_named_single && !.is_single_color_literal(palette)
    if (use_named_palette) {
      ggplot2::discrete_scale(
        aesthetics = "fill",
        palette = function(n) get_palette(name = palette, n = n, continuous = FALSE)
      )
    } else {
      pal <- as.character(palette)
      if (!is.null(names(palette)) && any(nzchar(names(palette)))) {
        ggplot2::scale_fill_manual(values = palette, drop = FALSE)
      } else {
        ggplot2::discrete_scale(
          aesthetics = "fill",
          palette = function(n) {
            if (length(pal) >= n) pal[seq_len(n)] else grDevices::colorRampPalette(pal)(n)
          }
        )
      }
    }
  }
}
