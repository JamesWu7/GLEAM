# GLEAM plotting theme

GLEAM plotting theme

## Usage

``` r
gleam_theme(
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
)
```

## Arguments

- base_size:

  Base text size used for the overall theme.

- title_size:

  Plot title size.

- subtitle_size:

  Plot subtitle size.

- axis_title_size:

  Axis title size.

- axis_text_size:

  Axis text size.

- legend_title_size:

  Legend title size.

- legend_text_size:

  Legend text size.

- strip_text_size:

  Facet strip text size.

- font_family:

  Font family.

- font_face:

  Font face.

- title_color:

  Title text color.

- subtitle_color:

  Subtitle text color.

- axis_text_color:

  Axis text color.

- axis_title_color:

  Axis title color.

- legend_text_color:

  Legend text color.

- legend_title_color:

  Legend title color.

- strip_text_color:

  Facet strip text color.

- text_color:

  Global text color fallback.

## Value

A ggplot2 theme object.
