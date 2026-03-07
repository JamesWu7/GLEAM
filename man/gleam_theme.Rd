% Manual Rd file
\name{gleam_theme}
\alias{gleam_theme}
\title{GLEAM plotting theme}
\usage{
gleam_theme(
  base_size = 13,
  title_size = base_size + 2,
  subtitle_size = base_size,
  axis_title_size = base_size,
  axis_text_size = max(8, base_size - 1),
  legend_title_size = base_size,
  legend_text_size = max(8, base_size - 1),
  strip_text_size = base_size,
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
}
\arguments{
\item{base_size}{Base text size used for the overall theme.}
\item{title_size}{Plot title size.}
\item{subtitle_size}{Plot subtitle size.}
\item{axis_title_size}{Axis title size.}
\item{axis_text_size}{Axis text size.}
\item{legend_title_size}{Legend title size.}
\item{legend_text_size}{Legend text size.}
\item{strip_text_size}{Facet strip text size.}
\item{font_family}{Font family.}
\item{font_face}{Font face.}
\item{title_color}{Title text color.}
\item{subtitle_color}{Subtitle text color.}
\item{axis_text_color}{Axis text color.}
\item{axis_title_color}{Axis title color.}
\item{legend_text_color}{Legend text color.}
\item{legend_title_color}{Legend title color.}
\item{strip_text_color}{Facet strip text color.}
\item{text_color}{Global text color fallback.}
}
\value{
A ggplot2 theme object.
}
\description{
Theme helper for readable, publication-friendly defaults with centralized
font, text-size, and color controls across GLEAM plots.
}
