#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-pipeline_outputs/ci_explorer_agent}"
mkdir -p "$OUT_DIR"

BUILD_LOG="$OUT_DIR/build.log"
CHECK_LOG="$OUT_DIR/check.log"

set +e
R CMD build . > "$BUILD_LOG" 2>&1
build_exit=$?
pkg_tar=$(ls -t ./*.tar.gz 2>/dev/null | head -n 1)
check_exit=99
if [[ $build_exit -eq 0 && -n "${pkg_tar:-}" ]]; then
  _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --as-cran "$pkg_tar" > "$CHECK_LOG" 2>&1
  check_exit=$?
else
  echo "Build failed; check step skipped." > "$CHECK_LOG"
fi
set -e

status_line=$(grep -E '^Status:' "$CHECK_LOG" | tail -n 1 || true)
error_line=$(grep -E '(^\* ERROR| ERROR$|Error: )' "$CHECK_LOG" | tail -n 1 || true)

{
  echo "build_exit=$build_exit"
  echo "check_exit=$check_exit"
} > "$OUT_DIR/status.env"

echo "$build_exit" > "$OUT_DIR/build_exit_code.txt"
echo "$check_exit" > "$OUT_DIR/check_exit_code.txt"

{
  echo "# ci_explorer_agent"
  echo
  echo "## Execution"
  echo "- build exit code: $build_exit"
  echo "- check exit code: $check_exit"
  echo
  if [[ -n "$status_line" ]]; then
    echo "## R CMD check status"
    echo "- $status_line"
    echo
  fi
  if [[ -n "$error_line" ]]; then
    echo "## Last error-like line"
    echo "- $error_line"
    echo
  fi
  echo "## Logs"
  echo "- build: \\`$BUILD_LOG\\`"
  echo "- check: \\`$CHECK_LOG\\`"
} > "$OUT_DIR/findings.md"

# Never fail this stage; downstream bug analyzer decides severity.
exit 0
