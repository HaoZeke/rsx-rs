#!/usr/bin/env bash
set -euo pipefail

ROOT=$(git rev-parse --show-toplevel)
SCRIPT="${RSX_PACKAGE_SCRIPT:-${ROOT}/repro/benchmarks/package_for_figshare.sh}"
OUTDIR=$(mktemp -d "${TMPDIR:-/tmp}/rsx-package-test.XXXXXX")

RSX_REPRO_OUTDIR="${OUTDIR}" bash "${SCRIPT}" >"${OUTDIR}/package.log"

shopt -s nullglob
packages=("${OUTDIR}"/rsx-benchmark-package-*)
shopt -u nullglob

if [[ "${#packages[@]}" -ne 1 ]]; then
    echo "expected exactly one package directory in ${OUTDIR}" >&2
    exit 1
fi

PACKAGE="${packages[0]}"

required_files=(
    "${PACKAGE}/pixi.toml"
    "${PACKAGE}/MANIFEST.txt"
    "${PACKAGE}/benchmarks/generate_data.py"
    "${PACKAGE}/benchmarks/run_benchmarks.sh"
    "${PACKAGE}/benchmarks/plot_benchmarks.py"
    "${PACKAGE}/benchmarks/results/benchmark_results.csv"
    "${PACKAGE}/benchmarks.org"
)

required_dirs=(
    "${PACKAGE}/benchmarks/data"
    "${PACKAGE}/benchmarks/results"
)

for path in "${required_files[@]}"; do
    if [[ ! -f "${path}" ]]; then
        echo "missing required file: ${path}" >&2
        exit 1
    fi
done

for path in "${required_dirs[@]}"; do
    if [[ ! -d "${path}" ]]; then
        echo "missing required directory: ${path}" >&2
        exit 1
    fi
done

if [[ -e "${PACKAGE}/data" || -e "${PACKAGE}/benchmark_results.csv" ]]; then
    echo "package must preserve repository-relative benchmarks/ paths" >&2
    exit 1
fi

grep -q "pixi.toml" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/data/" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/results/benchmark_results.csv" "${PACKAGE}/MANIFEST.txt"

(
    cd "${PACKAGE}"
    python3 benchmarks/plot_benchmarks.py \
        --input benchmarks/results/benchmark_results.csv \
        --output "${OUTDIR}/figures" >"${OUTDIR}/plot.log"
)

test -f "${OUTDIR}/figures/benchmark_summary.txt"
test -f "${OUTDIR}/figures/benchmark_speedups.svg"

echo "Package layout validated in ${PACKAGE}"
