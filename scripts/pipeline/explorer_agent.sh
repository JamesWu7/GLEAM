#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-pipeline_outputs/explorer_agent}"
mkdir -p "$OUT_DIR"

r_files=$(find R -type f -name '*.R' | wc -l | tr -d ' ')
rd_files=$(find man -type f -name '*.Rd' | wc -l | tr -d ' ')
test_files=$(find tests/testthat -type f -name '*.R' 2>/dev/null | wc -l | tr -d ' ')
workflow_files=$(find .github/workflows -type f -name '*.y*ml' | wc -l | tr -d ' ')

pathway_calls=$(rg -n "plot_(embedding_score|spatial_score|dot_bar|violin|split_violin|ridge|pseudotime_score|trajectory_score|box|pseudobulk_box)\\([^\\)]*pathway\\s*=" scripts vignettes README.md tests/testthat R 2>/dev/null || true)
old_methods=$(rg -n "zmean|singscore_like|ssgsea_like|method\\s*=\\s*\"auc\"|method\\s*=\\s*\"robust\"" R tests vignettes README.md man 2>/dev/null || true)
direct_push_main=$(rg -n "git push origin HEAD:\$\{GITHUB_REF_NAME\}" .github/workflows 2>/dev/null || true)

{
  echo "# explorer_agent"
  echo
  echo "## Repository snapshot"
  echo "- R files: $r_files"
  echo "- Rd files: $rd_files"
  echo "- testthat files: $test_files"
  echo "- workflow files: $workflow_files"
  echo
  echo "## Potential issues"

  if [[ -n "$pathway_calls" ]]; then
    echo "- [HIGH] Deprecated plot argument usage detected (pathway=)."
    echo
    echo '```text'
    echo "$pathway_calls"
    echo '```'
  else
    echo "- [OK] No deprecated pathway= plot argument usage detected."
  fi

  if [[ -n "$old_methods" ]]; then
    echo "- [MEDIUM] Legacy scoring method tokens detected."
    echo
    echo '```text'
    echo "$old_methods"
    echo '```'
  else
    echo "- [OK] No legacy scoring method tokens detected."
  fi

  if [[ -n "$direct_push_main" ]]; then
    echo "- [MEDIUM] Workflow step directly pushes to branch head (review for no-direct-push policy)."
    echo
    echo '```text'
    echo "$direct_push_main"
    echo '```'
  else
    echo "- [OK] No direct push-to-branch step detected in workflows."
  fi
} > "$OUT_DIR/findings.md"

printf 'status=ok\n' > "$OUT_DIR/status.env"
