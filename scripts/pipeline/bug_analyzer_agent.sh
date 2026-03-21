#!/usr/bin/env bash
set -euo pipefail

IN_DIR="${1:-pipeline_inputs}"
OUT_DIR="${2:-pipeline_outputs/bug_analyzer_agent}"
mkdir -p "$OUT_DIR"

ci_log="$IN_DIR/ci_explorer_agent/check.log"
pkgdown_log="$IN_DIR/pkgdown_explorer_agent/pkgdown.log"
explorer_md="$IN_DIR/explorer_agent/findings.md"
reviewer_md="$IN_DIR/reviewer_agent/findings.md"

combined_tmp="$OUT_DIR/combined.log"
{
  [[ -f "$ci_log" ]] && cat "$ci_log"
  [[ -f "$pkgdown_log" ]] && cat "$pkgdown_log"
  [[ -f "$explorer_md" ]] && cat "$explorer_md"
  [[ -f "$reviewer_md" ]] && cat "$reviewer_md"
} > "$combined_tmp"

error_type="UNKNOWN"
file_location="multiple"
root_cause="No known signature matched."
fix_strategy="Review logs manually and implement targeted patch."

if rg -n "unused argument \(pathway\s*=|unused argument \(pathway =" "$combined_tmp" >/dev/null 2>&1; then
  error_type="API_ARG_MISMATCH"
  file_location="scripts/*.Rmd, scripts/*.R"
  root_cause="Plot API migrated from pathway to signature, but call sites still pass pathway=."
  fix_strategy="Replace pathway= with signature= for plot_* calls and rerender homepage/tutorial artifacts."
elif rg -n "No geneset overlaps with expression matrix genes" "$combined_tmp" >/dev/null 2>&1; then
  error_type="SIGNED_GENESET_PARSE"
  file_location="R/geneset_utils.R"
  root_cause="Signed geneset list(up/down) was flattened incorrectly during coercion."
  fix_strategy="Preserve up/down structure in as_geneset(); add regression test."
elif rg -n "scale_name argument of .*discrete_scale\(\) is deprecated" "$combined_tmp" >/dev/null 2>&1; then
  error_type="GGPLOT_DEPRECATION"
  file_location="R/palette_utils.R"
  root_cause="Deprecated ggplot2 parameter scale_name still in use."
  fix_strategy="Remove scale_name from discrete_scale helper calls."
elif rg -n "Spatial metadata must contain coordinate columns|Failed to extract spatial coordinates" "$combined_tmp" >/dev/null 2>&1; then
  error_type="SPATIAL_COORD_FAILURE"
  file_location="R/spatial_utils.R, vignettes/GLEAM_full_spatial_workflow.Rmd"
  root_cause="Spatial coordinate extraction assumptions do not match provided object/meta structure."
  fix_strategy="Harden extract_spatial_coords() fallbacks and preflight-check in vignette before plotting."
elif rg -n "vignette builder 'knitr' not found" "$combined_tmp" >/dev/null 2>&1; then
  error_type="MISSING_VIGNETTE_BUILDER"
  file_location=".github/workflows/R-CMD-check.yaml"
  root_cause="Vignette build deps are missing in CI environment."
  fix_strategy="Install knitr/rmarkdown in check workflow before check-r-package step."
fi

{
  printf 'error_type=%q\n' "$error_type"
  printf 'file_location=%q\n' "$file_location"
  printf 'root_cause=%q\n' "$root_cause"
  printf 'fix_strategy=%q\n' "$fix_strategy"
} > "$OUT_DIR/bug_report.env"

{
  echo "# bug_analyzer_agent"
  echo
  echo "## Root cause summary"
  echo "- error_type: $error_type"
  echo "- file_location: $file_location"
  echo "- root_cause: $root_cause"
  echo "- fix_strategy: $fix_strategy"
} > "$OUT_DIR/bug_report.md"

exit 0
