#!/usr/bin/env bash
set -euo pipefail

IN_DIR="${1:-pipeline_inputs}"
OUT_DIR="${2:-pipeline_outputs/final_reviewer_agent}"
mkdir -p "$OUT_DIR"

bug_env="$IN_DIR/bug_analyzer_agent/bug_report.env"
validation_env="$IN_DIR/test_runner_agent/validation.env"
changes_md="$IN_DIR/fix_engineer_agent/changes.md"

error_type="UNKNOWN"
file_location="multiple"
root_cause="n/a"
fix_strategy="n/a"
overall="unknown"
build_exit="NA"
check_exit="NA"
pkgdown_exit="NA"

if [[ -f "$bug_env" ]]; then
  # shellcheck disable=SC1090
  source "$bug_env"
fi
if [[ -f "$validation_env" ]]; then
  # shellcheck disable=SC1090
  source "$validation_env"
fi

readiness="needs_changes"
if [[ "${overall:-unknown}" == "passed" ]]; then
  readiness="ready_for_human_review"
fi

{
  echo "readiness=$readiness"
  echo "overall=${overall:-unknown}"
  echo "error_type=${error_type:-UNKNOWN}"
} > "$OUT_DIR/review.env"

{
  echo "# final_reviewer_agent"
  echo
  echo "## Decision"
  echo "- readiness: $readiness"
  echo
  echo "## Regression check"
  echo "- build_exit: ${build_exit:-NA}"
  echo "- check_exit: ${check_exit:-NA}"
  echo "- pkgdown_exit: ${pkgdown_exit:-NA}"
  echo
  echo "## Root cause linkage"
  echo "- error_type: ${error_type:-UNKNOWN}"
  echo "- file_location: ${file_location:-multiple}"
  echo "- root_cause: ${root_cause:-n/a}"
  echo "- fix_strategy: ${fix_strategy:-n/a}"
  echo
  if [[ -f "$changes_md" ]]; then
    echo "## Applied patch summary"
    cat "$changes_md"
  fi
} > "$OUT_DIR/review.md"

exit 0
