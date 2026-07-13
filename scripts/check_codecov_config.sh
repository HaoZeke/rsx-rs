#!/usr/bin/env bash
# Structural gate: multi-flag Codecov wiring for rust/python/r stays intact.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
COV_YML="$ROOT/codecov.yml"
WF="$ROOT/.github/workflows/coverage.yml"
fail=0

die() { echo "ERROR: $*" >&2; fail=1; }
ok() { echo "OK: $*"; }

[[ -f "$COV_YML" ]] || die "missing $COV_YML"
[[ -f "$WF" ]] || die "missing $WF"

for flag in rust python r; do
  if grep -qE "name:[[:space:]]*${flag}" "$COV_YML"; then
    ok "codecov.yml flag $flag"
  else
    die "codecov.yml missing flag $flag"
  fi
done

grep -q 'informational: true' "$COV_YML" || die "codecov.yml missing informational: true"
grep -q 'carryforward: true' "$COV_YML" || die "codecov.yml missing carryforward: true"
ok "codecov.yml statuses/carryforward"

for flag in rust python r; do
  if grep -E "flags:[[:space:]]*${flag}" "$WF" | grep -vq '^\s*#'; then
    ok "coverage.yml upload flags: $flag"
  else
    die "coverage.yml missing active flags: $flag"
  fi
done

grep -q 'codecov/codecov-action' "$WF" || die "coverage.yml missing codecov-action"
grep -q 'fail_ci_if_error: false' "$WF" || die "coverage.yml missing fail_ci_if_error: false"
grep -q 'use_oidc: true' "$WF" || die "coverage.yml missing use_oidc: true"
grep -q 'id-token: write' "$WF" || die "coverage.yml missing id-token: write for OIDC"
grep -q 'app.codecov.io' "$WF" || die "coverage.yml missing app.codecov.io note"
ok "coverage.yml OIDC soft-fail + docs"

grep -qE 'cargo llvm-cov|run_coverage_rust\.sh' "$WF" || die "missing rust coverage generator"
grep -qE 'run_coverage_python\.sh|maturin' "$WF" || die "missing python coverage generator"
grep -qE 'run_coverage_r\.sh|covr' "$WF" || die "missing R coverage generator"
ok "real coverage generators referenced"

for s in run_coverage_rust.sh run_coverage_python.sh run_coverage_r.sh; do
  [[ -x "$ROOT/scripts/$s" ]] || die "missing executable scripts/$s"
  ok "scripts/$s"
done

if [[ "$fail" -ne 0 ]]; then
  echo "check_codecov_config: FAILED" >&2
  exit 1
fi
echo "check_codecov_config: all checks passed"
