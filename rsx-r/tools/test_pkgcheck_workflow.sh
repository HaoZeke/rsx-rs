#!/bin/sh
# Assert the committed rsxr-pkgcheck workflow meets the monorepo pkgcheck contract.
set -eu
ROOT="$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)"
WF="$ROOT/.github/workflows/rsxr-pkgcheck.yaml"
test -f "$WF"
grep -q 'name: rsxr-pkgcheck' "$WF"
grep -q 'stage_rsxcore.sh' "$WF"
grep -q 'dtolnay/rust-toolchain' "$WF"
grep -q 'zlib1g-dev' "$WF"
grep -q 'libbz2-dev' "$WF"
grep -q 'liblzma-dev' "$WF"
grep -q 'working-directory: rsx-r' "$WF"
grep -q 'pkgcheck::pkgcheck' "$WF"
grep -q 'workflow_dispatch' "$WF"
grep -q 'rsx-r/\*\*' "$WF"
grep -q 'rsxcore/\*\*' "$WF"
# Must not be Docker-only stock action without custom setup (filename may mention action in comments)
! grep -E 'uses:.*pkgcheck-action' "$WF"
grep -q 'rsxr-pkgcheck.yaml' "$ROOT/rsx-r/README.md"
grep -q 'Config/Needs/pkgcheck' "$ROOT/rsx-r/DESCRIPTION"
echo "test_pkgcheck_workflow.sh: all assertions passed"
