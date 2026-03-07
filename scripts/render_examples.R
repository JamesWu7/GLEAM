out_dir <- file.path("docs", "tutorial_html")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  message("[GLEAM] Package 'rmarkdown' is not installed; skip local HTML rendering.")
} else {
  vigs <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
  if (length(vigs) == 0) {
    message("[GLEAM] No vignette Rmd files found.")
  } else {
    preferred <- c(
      "GLEAM_quickstart_scrna.Rmd",
      "GLEAM_quickstart_spatial.Rmd",
      "GLEAM_custom_geneset.Rmd",
      "GLEAM_supported_genesets.Rmd",
      "GLEAM_visualization_parameters.Rmd",
      "GLEAM_function_categories.Rmd",
      "GLEAM_full_scrna_workflow.Rmd",
      "GLEAM_full_spatial_workflow.Rmd",
      "GLEAM_citation.Rmd"
    )
    ord <- match(basename(vigs), preferred)
    vigs <- c(vigs[order(is.na(ord), ord)])

    message("[GLEAM] Rendering ", length(vigs), " tutorial files to ", out_dir)
    for (vf in vigs) {
      message("[GLEAM] Rendering: ", vf)
      tryCatch(
        rmarkdown::render(vf, output_format = "html_document", output_dir = out_dir, quiet = TRUE),
        error = function(e) message("[GLEAM] Skipped (error): ", e$message)
      )
    }
  }
}
