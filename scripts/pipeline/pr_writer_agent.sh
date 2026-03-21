#!/usr/bin/env bash
set -euo pipefail

IN_DIR="${1:-pipeline_inputs}"
OUT_DIR="${2:-pipeline_outputs/pr_writer_agent}"
mkdir -p "$OUT_DIR"

bug_env="$IN_DIR/bug_analyzer_agent/bug_report.env"
validation_env="$IN_DIR/test_runner_agent/validation.env"
review_env="$IN_DIR/final_reviewer_agent/review.env"

error_type="UNKNOWN"
file_location="multiple"
root_cause="n/a"
fix_strategy="n/a"
overall="unknown"
readiness="needs_changes"

if [[ -f "$bug_env" ]]; then
  # shellcheck disable=SC1090
  source "$bug_env"
fi
if [[ -f "$validation_env" ]]; then
  # shellcheck disable=SC1090
  source "$validation_env"
fi
if [[ -f "$review_env" ]]; then
  # shellcheck disable=SC1090
  source "$review_env"
fi

case "${error_type:-UNKNOWN}" in
  API_ARG_MISMATCH)
    pr_title="fix(pkgdown): migrate plot args from pathway to signature"
    commit_msg="fix(pkgdown): align plot argument names with signature API"
    ;;
  SIGNED_GENESET_PARSE)
    pr_title="fix(scoring): preserve signed up/down geneset structure"
    commit_msg="fix(scoring): keep signed genesets intact in as_geneset"
    ;;
  GGPLOT_DEPRECATION)
    pr_title="chore(plot): remove deprecated ggplot2 discrete_scale arg"
    commit_msg="chore(plot): drop deprecated scale_name argument"
    ;;
  SPATIAL_COORD_FAILURE)
    pr_title="fix(spatial): harden coordinate extraction fallbacks"
    commit_msg="fix(spatial): add robust spatial coordinate fallback paths"
    ;;
  *)
    pr_title="chore(ci): autonomous maintenance pipeline updates"
    commit_msg="chore(ci): apply autonomous pipeline fixes"
    ;;
esac

cat > "$OUT_DIR/pr_title.txt" <<TXT
$pr_title
TXT

cat > "$OUT_DIR/commit_message.txt" <<TXT
$commit_msg
TXT

cat > "$OUT_DIR/pr_body.md" <<MD
## Problem
Autonomous pipeline detected issues during analysis/validation and prepared a bounded fix patch for human review.

## Root Cause
- **error_type**: ${error_type:-UNKNOWN}
- **file_location**: ${file_location:-multiple}
- **root_cause**: ${root_cause:-n/a}

## Solution
- Applied minimal, module-bounded changes using deterministic fix rules.
- Preserved API compatibility where possible and avoided direct push to \\`main\\`.

## Validation
- **overall**: ${overall:-unknown}
- Validation artifacts are attached to the workflow run (analysis, fix, validation, final review).

## User Impact
- Improves CI/pkgdown stability and reduces recurring maintenance failures.
- Keeps workflow review-driven: PR is created for GitHub human review before merge.

## Human Review Checklist (GitHub)
1. Open the PR Checks tab and inspect uploaded artifacts.
2. Verify pkgdown outputs and docs links in the preview/deployed site.
3. Verify key tutorials and plots (UMAP, spatial, dot-bar, violin, ridge).
4. Confirm API consistency and no regressions.
5. Merge only after explicit approval.

## Pipeline State
- **readiness**: ${readiness:-needs_changes}
MD

exit 0
