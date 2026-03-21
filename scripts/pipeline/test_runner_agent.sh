#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:-pipeline_outputs/test_runner_agent}"
mkdir -p "$OUT_DIR"

CHECK_LOG="$OUT_DIR/check.log"
PKGDOWN_LOG="$OUT_DIR/pkgdown.log"

set +e
R CMD build . > "$OUT_DIR/build.log" 2>&1
build_exit=$?
pkg_tar=$(ls -t ./*.tar.gz 2>/dev/null | head -n 1)
check_exit=99
if [[ $build_exit -eq 0 && -n "${pkg_tar:-}" ]]; then
  _R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --as-cran "$pkg_tar" > "$CHECK_LOG" 2>&1
  check_exit=$?
else
  echo "Build failed; check skipped." > "$CHECK_LOG"
fi

scripts/pipeline/pkgdown_build_and_validate.sh > "$PKGDOWN_LOG" 2>&1
pkgdown_exit=$?
set -e

overall="passed"
if [[ $build_exit -ne 0 || $check_exit -ne 0 || $pkgdown_exit -ne 0 ]]; then
  overall="failed"
fi

{
  echo "build_exit=$build_exit"
  echo "check_exit=$check_exit"
  echo "pkgdown_exit=$pkgdown_exit"
  echo "overall=$overall"
} > "$OUT_DIR/validation.env"

{
  echo "# test_runner_agent"
  echo
  echo "## Validation summary"
  echo "- build_exit: $build_exit"
  echo "- check_exit: $check_exit"
  echo "- pkgdown_exit: $pkgdown_exit"
  echo "- overall: $overall"
  echo
  echo "## Notes"
  echo "- This stage is non-blocking for PR generation; failures are surfaced in PR body for human review."
} > "$OUT_DIR/summary.md"

# Keep stage non-blocking; final reviewer decides readiness.
exit 0
