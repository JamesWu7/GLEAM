#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-pipeline_outputs/reviewer_agent}"
mkdir -p "$OUT_DIR"

man_pathway=$(rg -n "\\\item\\{pathway\\}" man/plot_*.Rd 2>/dev/null || true)
roxy_stale=$(rg -n "Please edit documentation in R/score_zmean\\.R" man 2>/dev/null || true)
legacy_alias=$(rg -n "Legacy alias|legacy argument name" man R 2>/dev/null || true)

{
  echo "# reviewer_agent"
  echo
  echo "## API / docs consistency review"
  if [[ -n "$man_pathway" ]]; then
    echo "- [MEDIUM] Plot Rd files still exposing old `pathway` argument."
    echo
    echo '```text'
    echo "$man_pathway"
    echo '```'
  else
    echo "- [OK] Plot Rd argument names are aligned."
  fi

  if [[ -n "$roxy_stale" ]]; then
    echo "- [LOW] Roxygen source hint references removed file(s)."
    echo
    echo '```text'
    echo "$roxy_stale"
    echo '```'
  else
    echo "- [OK] No stale roxygen source hints found."
  fi

  if [[ -n "$legacy_alias" ]]; then
    echo "- [LOW] Residual legacy-compatibility wording found in code/docs."
    echo
    echo '```text'
    echo "$legacy_alias"
    echo '```'
  else
    echo "- [OK] No residual legacy-compatibility wording found."
  fi
} > "$OUT_DIR/findings.md"

printf 'status=ok\n' > "$OUT_DIR/status.env"
