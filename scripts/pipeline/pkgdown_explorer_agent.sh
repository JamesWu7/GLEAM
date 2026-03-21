#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-pipeline_outputs/pkgdown_explorer_agent}"
mkdir -p "$OUT_DIR"

LOG_FILE="$OUT_DIR/pkgdown.log"

set +e
scripts/pipeline/pkgdown_build_and_validate.sh > "$LOG_FILE" 2>&1
pkgdown_exit=$?
set -e

{
  echo "pkgdown_exit=$pkgdown_exit"
  if [[ $pkgdown_exit -eq 0 ]]; then
    echo "docs_complete=true"
  else
    echo "docs_complete=false"
  fi
} > "$OUT_DIR/status.env"

echo "$pkgdown_exit" > "$OUT_DIR/pkgdown_exit_code.txt"

last_error=$(grep -E 'Error:|Failed to render|Execution halted|unused argument' "$LOG_FILE" | tail -n 1 || true)

{
  echo "# pkgdown_explorer_agent"
  echo
  echo "## Execution"
  echo "- pkgdown exit code: $pkgdown_exit"
  if [[ $pkgdown_exit -eq 0 ]]; then
    echo "- docs key pages: present"
  else
    echo "- docs key pages: missing or pkgdown validation failed"
  fi
  echo
  if [[ -n "$last_error" ]]; then
    echo "## Last error-like line"
    echo "- $last_error"
    echo
  fi
  echo "## Log"
  echo "- \\`$LOG_FILE\\`"
} > "$OUT_DIR/findings.md"

exit 0
