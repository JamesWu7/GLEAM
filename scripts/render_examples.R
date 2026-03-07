out_dir <- file.path("docs", "tutorial_html")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vigs <- list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE)
if (length(vigs) == 0) {
  message("[GLEAM] No vignette Rmd files found.")
} else {
  message("[GLEAM] Rendering ", length(vigs), " tutorial files to ", out_dir)
  for (vf in vigs) {
    message("[GLEAM] Rendering: ", vf)
    tryCatch(
      rmarkdown::render(vf, output_format = "html_document", output_dir = out_dir, quiet = TRUE),
      error = function(e) message("[GLEAM] Skipped (error): ", e$message)
    )
  }
}
