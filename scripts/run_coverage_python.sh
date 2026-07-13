#!/usr/bin/env bash
# Instrument pyrsx PyO3 (rsx-python/src) via cargo-llvm-cov + maturin + pytest,
# and collect pure-Python coverage under rsx-python/python/.
#
# Do NOT clear RUSTC_WRAPPER after show-env — that strips instrumentation and
# yields zero .profraw (same bug class as readcon-core).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT_LCOV="${1:-python_lcov.info}"
OUT_XML="${2:-python_coverage.xml}"
VENV="${RSX_PY_VENV:-$ROOT/.venv-coverage-python}"

unset SCCACHE_GHA_ENABLED || true
if [[ "${RUSTC_WRAPPER:-}" == *sccache* ]]; then
  unset RUSTC_WRAPPER || true
fi
export CARGO_INCREMENTAL=0

if [[ -z "${CC:-}" ]] && command -v clang >/dev/null 2>&1; then
  export CC=clang CXX="${CXX:-clang++}"
fi
if [[ -z "${PYO3_PYTHON:-}" ]] && command -v python3 >/dev/null 2>&1; then
  export PYO3_PYTHON="$(command -v python3)"
fi

echo "==> python binding coverage (pyrsx)"
cargo llvm-cov clean --workspace 2>/dev/null || true

# shellcheck disable=SC1090
source <(cargo llvm-cov show-env --sh 2>/dev/null || cargo llvm-cov show-env --export-prefix)

if [[ -z "${RUSTC_WRAPPER:-}" ]] || [[ "${RUSTC_WRAPPER}" == *sccache* ]]; then
  echo "ERROR: RUSTC_WRAPPER must be cargo-llvm-cov after show-env (got: ${RUSTC_WRAPPER:-empty})" >&2
  exit 1
fi
echo "    RUSTC_WRAPPER=${RUSTC_WRAPPER}"
echo "    LLVM_PROFILE_FILE=${LLVM_PROFILE_FILE:-unset}"

if [[ ! -d "$VENV" ]]; then
  "${PYO3_PYTHON:-python3}" -m venv "$VENV"
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install -U pip -q
# pandas: test_marker_table_from_arrow; polars: backend-agnostic high-level tests
python -m pip install maturin pytest pytest-cov coverage \
  click 'narwhals>=1.0' 'pyarrow>=14' pandas polars -q

case " ${RUSTFLAGS:-} " in
  *" -C link-arg=-fuse-ld=bfd "*) ;;
  *) export RUSTFLAGS="${RUSTFLAGS:-} -C link-arg=-fuse-ld=bfd" ;;
esac

(
  cd "$ROOT/rsx-python"
  maturin develop
)

python -m pytest "$ROOT/rsx-python/tests" -q --tb=short \
  --cov=pyrsx \
  --cov-report=xml:"$ROOT/$OUT_XML" \
  --cov-report=term-missing \
  --cov-branch

PROF_GLOB="${CARGO_LLVM_COV_TARGET_DIR:-$ROOT/target}"
if ! compgen -G "${PROF_GLOB}"/*.profraw > /dev/null 2>&1 \
   && ! compgen -G "${ROOT}/target"/*.profraw > /dev/null 2>&1; then
  echo "    no profraw yet; probing import flush" >&2
  python -c 'import pyrsx; print("pyrsx", getattr(pyrsx, "__file__", "?"))'
fi

TMP_LCOV="$(mktemp)"
cargo llvm-cov report --lcov --output-path "$TMP_LCOV"
python3 - "$TMP_LCOV" "$ROOT/$OUT_LCOV" <<'PY'
import sys
inp, outp = sys.argv[1], sys.argv[2]
keep = False
buf = []
out = []
for line in open(inp, encoding="utf-8", errors="replace"):
    if line.startswith("SF:"):
        if buf and keep:
            out.extend(buf)
        buf = [line]
        path = line[3:].strip().replace("\\", "/")
        keep = (
            "/rsx-python/src/" in path
            or path.endswith("rsx-python/src/lib.rs")
            or ("/rsx-python/" in path and path.endswith(".rs"))
        )
    elif line.startswith("end_of_record"):
        buf.append(line)
        if keep:
            out.extend(buf)
        buf = []
        keep = False
    else:
        buf.append(line)
if buf and keep:
    out.extend(buf)
text = "".join(out)
if "SF:" not in text:
    raise SystemExit(
        "no rsx-python/src records in llvm-cov report — was pyrsx instrumented?"
    )
open(outp, "w", encoding="utf-8").write(text)
hits = tot = 0
for line in text.splitlines():
    if line.startswith("DA:"):
        h = int(line.split(":")[1].split(",")[1])
        tot += 1
        if h > 0:
            hits += 1
print(f"pyrsx rust lcov {100 * hits / tot:.1f}% {hits}/{tot}")
PY
rm -f "$TMP_LCOV"
test -s "$ROOT/$OUT_LCOV"
test -s "$ROOT/$OUT_XML"
echo "OK pure-python coverage XML: $ROOT/$OUT_XML"
echo "OK wrote $ROOT/$OUT_LCOV"
