#!/usr/bin/env bash
set -euo pipefail

IN_DIR="${1:-pipeline_inputs/bug_analyzer_agent}"
OUT_DIR="${2:-pipeline_outputs/fix_engineer_agent}"
mkdir -p "$OUT_DIR"

bug_env="$IN_DIR/bug_report.env"
if [[ -f "$bug_env" ]]; then
  # shellcheck disable=SC1090
  source "$bug_env"
fi

# Rule 1: migrate deprecated plot arg pathway= -> signature= at known call sites.
mapfile -t target_files < <(find scripts vignettes tests README.md README.Rmd -type f 2>/dev/null | sed '/\/\.git\//d')
if [[ ${#target_files[@]} -gt 0 ]]; then
  perl -0pi -e 's/(plot_(?:embedding_score|spatial_score|dot_bar|violin|split_violin|ridge|pseudotime_score|trajectory_score|box|pseudobulk_box)\s*\([^\)]*?)\bpathway\s*=/\1signature =/gs' "${target_files[@]}" || true
fi

# Rule 2: remove deprecated ggplot2 discrete_scale(scale_name=...) argument.
if [[ -f "R/palette_utils.R" ]]; then
  perl -0pi -e 's/,\s*\n\s*scale_name\s*=\s*"gleam_(?:color|fill)"\s*//g' R/palette_utils.R || true
fi

# Rule 3: enforce no direct push-to-main in pkgdown workflow by using Pages artifact deploy model.
if [[ -f ".github/workflows/pkgdown.yaml" ]]; then
  if rg -n "git push origin HEAD:\$\{GITHUB_REF_NAME\}" .github/workflows/pkgdown.yaml >/dev/null 2>&1; then
    perl -0pi -e 's/\n\s*- name: Commit docs to main\n\s*run:\s*\|\n(?:\s{10}.*\n)+/\n/gs' .github/workflows/pkgdown.yaml || true
  fi
fi

if git diff --quiet; then
  echo "changes_detected=false" > "$OUT_DIR/status.env"
  : > "$OUT_DIR/autofix.patch"
  {
    echo "# fix_engineer_agent"
    echo
    echo "No code changes were required for the detected issue profile."
  } > "$OUT_DIR/changes.md"
  exit 0
fi

# Keep patch as artifact for downstream validation and PR creation.
git diff --binary > "$OUT_DIR/autofix.patch"
changed_files=$(git diff --name-only)

{
  echo "changes_detected=true"
  if [[ -n "${error_type:-}" ]]; then
    echo "error_type=${error_type}"
  fi
} > "$OUT_DIR/status.env"

{
  echo "# fix_engineer_agent"
  echo
  echo "## Applied changes"
  echo '```text'
  echo "$changed_files"
  echo '```'
} > "$OUT_DIR/changes.md"

exit 0
