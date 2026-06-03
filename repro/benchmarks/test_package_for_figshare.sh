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
    "${PACKAGE}/benchmarks/literature_datasets.tsv"
    "${PACKAGE}/benchmarks/run_literature_benchmarks.py"
    "${PACKAGE}/benchmarks/plot_literature_benchmarks.py"
    "${PACKAGE}/benchmarks/plot_prior_sensitivity_heatmap.py"
    "${PACKAGE}/scripts/prepare_python_build.py"
    "${PACKAGE}/benchmarks/slurm/literature_biology_prior_sensitivity.sbatch"
    "${PACKAGE}/benchmarks/slurm/literature_biology_low_depth.sbatch"
    "${PACKAGE}/benchmarks/results/benchmark_results.csv"
    "${PACKAGE}/benchmarks/results/literature_benchmark_results.csv"
    "${PACKAGE}/benchmarks/results/literature_speed_comparison.csv"
    "${PACKAGE}/benchmarks/results/literature_binding_results.csv"
    "${PACKAGE}/benchmarks/results/literature_depth_stability.csv"
    "${PACKAGE}/benchmarks/results/prior_sensitivity_from_triage.csv"
    "${PACKAGE}/benchmarks/results/slurm/literature_speed_comparison_danio_albolineatus_manual.csv"
    "${PACKAGE}/benchmarks/results/slurm/triage_danio_albolineatus_pi0.001_psex0.8.tsv"
    "${PACKAGE}/docs/figures/literature_radsex_speedups.svg"
    "${PACKAGE}/docs/figures/literature_depth_stability.svg"
    "${PACKAGE}/benchmarks.org"
    "${PACKAGE}/literature_benchmarks.org"
)

required_dirs=(
    "${PACKAGE}/benchmarks/data"
    "${PACKAGE}/benchmarks/results"
    "${PACKAGE}/benchmarks/slurm"
    "${PACKAGE}/docs/figures"
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

if [[ -e "${PACKAGE}/benchmarks/literature-workdir" ]]; then
    echo "package must not include downloaded FASTQ workdir payloads" >&2
    exit 1
fi

grep -q "pixi.toml" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/data/" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/results/benchmark_results.csv" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/results/literature_speed_comparison.csv" "${PACKAGE}/MANIFEST.txt"
grep -q "benchmarks/results/slurm/triage_" "${PACKAGE}/MANIFEST.txt"
grep -q "docs/figures/literature_radsex_speedups.svg" "${PACKAGE}/MANIFEST.txt"
grep -q "scripts/prepare_python_build.py" "${PACKAGE}/MANIFEST.txt"

(
    cd "${PACKAGE}"
    python3 benchmarks/plot_benchmarks.py \
        --input benchmarks/results/benchmark_results.csv \
        --output "${OUTDIR}/figures" >"${OUTDIR}/plot.log"
)

test -f "${OUTDIR}/figures/benchmark_summary.txt"
test -f "${OUTDIR}/figures/benchmark_speedups.svg"

echo "Package layout validated in ${PACKAGE}"
