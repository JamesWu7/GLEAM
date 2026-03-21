#!/usr/bin/env bash
set -euo pipefail

Rscript -e 'pkgdown::clean_site(force = TRUE)'
Rscript -e 'source("scripts/generate_homepage_figures.R")'
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
