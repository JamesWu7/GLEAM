#!/usr/bin/env bash
set -euo pipefail

required_figs=(
  "assets/figures/embedding_signature_feature.png"
  "assets/figures/spatial_slice_signature.png"
  "assets/figures/signature_dotbar_compare.png"
)

Rscript -e 'pkgdown::clean_site(force = TRUE)'
if ! Rscript -e 'source("scripts/generate_homepage_figures.R")'; then
  echo "[GLEAM] WARNING: homepage figure generation failed; falling back to existing assets." >&2
fi

for f in "${required_figs[@]}"; do
  if [[ ! -f "${f}" ]]; then
    echo "[GLEAM] ERROR: required figure missing after generation fallback: ${f}" >&2
    exit 1
  fi
done

Rscript -e 'pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE)'

mkdir -p docs/assets
cp -R assets/. docs/assets/
mkdir -p docs/man/figures
cp -R man/figures/. docs/man/figures/

test -f docs/index.html
test -f docs/articles/index.html
test -f docs/articles/GLEAM_full_scrna_workflow.html
test -f docs/articles/GLEAM_full_spatial_workflow.html
test -f docs/articles/GLEAM_citation.html
test -f docs/reference/index.html
test -f docs/.nojekyll
