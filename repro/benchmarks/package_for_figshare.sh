#!/usr/bin/env bash
set -euo pipefail

# Package the benchmark data and results for DOI deposition.

ROOT=$(git rev-parse --show-toplevel)
OUTDIR="${RSX_REPRO_OUTDIR:-${ROOT}/repro/benchmarks/figshare_package}"
PACKAGE="${OUTDIR}/rsx-benchmark-package-$(git rev-parse --short HEAD)"

required_paths=(
    "${ROOT}/benchmarks/data"
    "${ROOT}/benchmarks/results/benchmark_results.csv"
    "${ROOT}/benchmarks/results/literature_speed_comparison.csv"
    "${ROOT}/benchmarks/results/prior_sensitivity_from_triage.csv"
    "${ROOT}/benchmarks/results/slurm/triage_danio_albolineatus_pi0.001_psex0.8.tsv"
    "${ROOT}/benchmarks/literature_datasets.tsv"
    "${ROOT}/benchmarks/slurm/literature_biology_prior_sensitivity.sbatch"
    "${ROOT}/docs/figures/literature_radsex_speedups.svg"
    "${ROOT}/pixi.toml"
    "${ROOT}/benchmarks/generate_data.py"
    "${ROOT}/benchmarks/run_benchmarks.sh"
    "${ROOT}/benchmarks/plot_benchmarks.py"
    "${ROOT}/repro/benchmarks.org"
    "${ROOT}/repro/literature_benchmarks.org"
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
mkdir -p "${PACKAGE}/docs"

echo "Copying benchmark inputs and results used in the manuscript..."
cp "${ROOT}/pixi.toml" "${PACKAGE}/"
cp -R "${ROOT}/benchmarks/data" "${PACKAGE}/benchmarks/"
cp -R "${ROOT}/benchmarks/results" "${PACKAGE}/benchmarks/"
cp -R "${ROOT}/benchmarks/slurm" "${PACKAGE}/benchmarks/"
cp -R "${ROOT}/docs/figures" "${PACKAGE}/docs/"
cp "${ROOT}/benchmarks/"*.py "${PACKAGE}/benchmarks/"
cp "${ROOT}/benchmarks/"*.sh "${PACKAGE}/benchmarks/"
cp "${ROOT}/benchmarks/literature_datasets.tsv" "${PACKAGE}/benchmarks/"
cp "${ROOT}/repro/benchmarks.org" "${PACKAGE}/"
cp "${ROOT}/repro/literature_benchmarks.org" "${PACKAGE}/"

# Create a manifest
cat > "${PACKAGE}/MANIFEST.txt" << EOF
rsx BMC Bioinformatics benchmark package
Generated on: $(date -u +%Y-%m-%dT%H:%M:%SZ)
Git commit: $(git rev-parse HEAD)
Paper: rsx software article (BMC Bioinformatics, 2026)

Contents:
- pixi.toml : Reproducible command environment definition
- benchmarks/data/ : Synthetic RAD-seq inputs and generated marker tables used for timings
- benchmarks/results/benchmark_results.csv : Synthetic benchmark CSV used by the regression figures
- benchmarks/results/literature_benchmark_results.csv : Published RADSex workflow dataset timings
- benchmarks/results/literature_speed_comparison.csv : Same-input rsx-rs vs C++ RADSex timings
- benchmarks/results/literature_binding_results.csv : Python binding feature timings on published marker tables
- benchmarks/results/literature_depth_stability.csv : Low-depth Bayesian and strict-call stability summary
- benchmarks/results/prior_sensitivity_from_triage.csv : Prior-sensitivity summary from the full triage grid
- benchmarks/results/slurm/triage_*_pi*_psex*.tsv : Per-cell Bayesian prior-sensitivity outputs
- benchmarks/results/slurm/literature_*_manual.csv : Per-dataset Slurm shards used by the aggregate tables
- docs/figures/literature_*.svg|pdf : Paper figures regenerated from the tracked result CSV/TSV files
- docs/figures/literature_radsex_speedups.svg : Same-input C++ RADSex vs rsx-rs speedup figure
- docs/figures/literature_depth_stability.svg : Low-depth strict/posterior/Bayes-factor stability figure
- benchmarks/*.py, benchmarks/*.sh, benchmarks/slurm/*.sbatch : Scripts to regenerate equivalent data
- benchmarks.org, literature_benchmarks.org : Full reproduction instructions

To reproduce the synthetic numbers, use the archived benchmarks/data/ and CSV directly,
or regenerate with the provided pixi environment and scripts (see benchmarks.org).
To reproduce the literature panel, use the public accessions in benchmarks/literature_datasets.tsv
and the Slurm/Pixi commands in literature_benchmarks.org.
EOF

echo "Package prepared in ${PACKAGE}"
echo "Upload the contents of ${PACKAGE} to a DOI repository as a single archive."
