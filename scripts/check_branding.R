targets <- c("README.md", "README.Rmd", "DESCRIPTION", "CONTRIBUTING.md", "DEVELOPMENT.md", "ROADMAP.md", "vignettes", ".github")
hits <- unlist(lapply(targets, function(p) {
  if (!file.exists(p)) return(character(0))
  suppressWarnings(system2("rg", c("-n", "-e", "scPathway", "-e", "scpathway", p), stdout = TRUE, stderr = FALSE))
}), use.names = FALSE)
hits <- unique(hits)
if (length(hits) == 0) {
  message("[GLEAM] Branding check passed (docs/repo presentation are GLEAM-only).")
} else {
  message("[GLEAM] Branding leftovers detected in presentation docs:")
  cat(paste(hits, collapse = "\n"), "\n")
}
