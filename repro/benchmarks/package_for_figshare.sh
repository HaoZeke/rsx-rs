#!/usr/bin/env bash
set -euo pipefail

# Package the benchmark data and results for DOI deposition.

ROOT=$(git rev-parse --show-toplevel)
OUTDIR="${RSX_REPRO_OUTDIR:-${ROOT}/repro/benchmarks/figshare_package}"
PACKAGE="${OUTDIR}/rsx-benchmark-package-$(git rev-parse --short HEAD)"

required_paths=(
    "${ROOT}/benchmarks/data"
    "${ROOT}/benchmarks/results/benchmark_results.csv"
    "${ROOT}/pixi.toml"
    "${ROOT}/benchmarks/generate_data.py"
    "${ROOT}/benchmarks/run_benchmarks.sh"
    "${ROOT}/benchmarks/plot_benchmarks.py"
    "${ROOT}/repro/benchmarks.org"
)

for path in "${required_paths[@]}"; do
    if [[ ! -e "${path}" ]]; then
        echo "missing required path: ${path}" >&2
        exit 1
    fi
done

if [[ -e "${PACKAGE}" ]]; then
    echo "package already exists: ${PACKAGE}" >&2
    exit 1
fi

mkdir -p "${PACKAGE}"
mkdir -p "${PACKAGE}/benchmarks"

echo "Copying benchmark inputs and results used in the manuscript..."
cp "${ROOT}/pixi.toml" "${PACKAGE}/"
cp -R "${ROOT}/benchmarks/data" "${PACKAGE}/benchmarks/"
mkdir -p "${PACKAGE}/benchmarks/results"
cp "${ROOT}/benchmarks/results/benchmark_results.csv" "${PACKAGE}/benchmarks/results/"
cp "${ROOT}/benchmarks/generate_data.py" "${PACKAGE}/benchmarks/"
cp "${ROOT}/benchmarks/run_benchmarks.sh" "${PACKAGE}/benchmarks/"
cp "${ROOT}/benchmarks/plot_benchmarks.py" "${PACKAGE}/benchmarks/"
cp "${ROOT}/repro/benchmarks.org" "${PACKAGE}/"

# Create a manifest
cat > "${PACKAGE}/MANIFEST.txt" << EOF
rsx BMC Bioinformatics benchmark package
Generated on: $(date -u +%Y-%m-%dT%H:%M:%SZ)
Git commit: $(git rev-parse HEAD)
Paper: rsx software article (BMC Bioinformatics, 2026)

Contents:
- pixi.toml : Reproducible command environment definition
- benchmarks/data/ : Synthetic RAD-seq inputs and generated marker tables used for timings
- benchmarks/results/benchmark_results.csv : The exact CSV used for all numbers and figures in the manuscript
- benchmarks/generate_data.py, benchmarks/run_benchmarks.sh, benchmarks/plot_benchmarks.py : Scripts to (re)generate equivalent data
- benchmarks.org : Full reproduction instructions

To reproduce the numbers in the paper, use the archived benchmarks/data/ and CSV directly,
or regenerate with the provided pixi environment and scripts (see benchmarks.org).
EOF

echo "Package prepared in ${PACKAGE}"
echo "Upload the contents of ${PACKAGE} to a DOI repository as a single archive."
