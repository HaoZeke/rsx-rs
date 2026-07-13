#!/usr/bin/env bash
# Produce Codecov LCOV (+ optional codecov JSON) for the Rust library surface.
#
# Primary package: rsxcore (includes c_api FFI). CLI binary is optional smoke.
# PyO3 (pyrsx) is covered under the python flag, not here.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT_LCOV="${1:-lcov.info}"
OUT_JSON="${2:-rust_codecov.json}"
# map pulls minimap2 (default in rsxcore); arrow* for high-level API paths.
# mpi intentionally omitted (needs OpenMPI + special runners).
FEATURES="${RSX_COV_FEATURES:-parallel,map,arrow-output,arrow-input,parquet-io}"
IGNORE='(/rsx-python/|/tests/|/benches/|/benchmarks/)'

unset RUSTC_WRAPPER SCCACHE_GHA_ENABLED || true
export RUSTC_WRAPPER=""
export CARGO_INCREMENTAL=0

if [[ -z "${CC:-}" ]] && command -v clang >/dev/null 2>&1; then
  export CC=clang CXX="${CXX:-clang++}"
fi

echo "==> cargo llvm-cov -p rsxcore (features=${FEATURES})"
# shellcheck disable=SC2086
cargo llvm-cov -p rsxcore --features ${FEATURES} \
  --no-fail-fast \
  --ignore-filename-regex="${IGNORE}" \
  --lcov --output-path "${OUT_LCOV}"

cargo llvm-cov report --codecov --output-path "${OUT_JSON}" \
  --ignore-filename-regex="${IGNORE}"

test -s "$OUT_LCOV"
test -s "$OUT_JSON"

python3 - "$OUT_LCOV" <<'PY'
import sys
path = sys.argv[1]
hits = tot = 0
for line in open(path, encoding="utf-8", errors="replace"):
    if line.startswith("DA:"):
        h = int(line.strip().split(":")[1].split(",")[1])
        tot += 1
        if h > 0:
            hits += 1
pct = 100 * hits / tot if tot else 0.0
print(f"rust lcov line coverage {pct:.2f}%  {hits}/{tot}  ({path})")
if tot == 0:
    sys.exit("empty LCOV")
PY
echo "OK wrote $OUT_LCOV and $OUT_JSON"
