#!/usr/bin/env bash
# literature_biology_prior_sensitivity.sbatch
# Re-run rsx triage on real literature marker tables with varied --prior and --linked-prob.
# Produces the data for prior sensitivity analysis (vissue internal-issue).
#
# Usage on elja (example):
#   RSX_DATASETS="danio_albolineatus notothenia_rossii plecoglossus_altivelis tinca_tinca" \
#   RSX_PRIORS="0.001 0.005 0.01 0.02 0.05" \
#   RSX_PSEX="0.8 0.85 0.9 0.95" \
#   sbatch --array=0-3 benchmarks/slurm/literature_biology_prior_sensitivity.sbatch
#
#SBATCH --job-name=rsx-bio-prior
#SBATCH --output=rsx-bio-prior-%A_%a.out
#SBATCH --error=rsx-bio-prior-%A_%a.err
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --export=ALL

set -euo pipefail

repo=${RSX_REPO:-$PWD}
cd "$repo"

if [[ -n "${PIXI_BIN:-}" ]]; then
    pixi=$PIXI_BIN
elif command -v pixi >/dev/null 2>&1; then
    pixi=$(command -v pixi)
elif [[ -x "$HOME/.pixi/bin/pixi" ]]; then
    pixi="$HOME/.pixi/bin/pixi"
else
    echo "pixi not found; set PIXI_BIN" >&2
    exit 127
fi

read -r -a datasets <<< "${RSX_DATASETS:-danio_albolineatus notothenia_rossii plecoglossus_altivelis tinca_tinca}"
index=${SLURM_ARRAY_TASK_ID:-0}
if (( index < 0 || index >= ${#datasets[@]} )); then
    echo "SLURM_ARRAY_TASK_ID=$index is outside ${#datasets[@]} datasets" >&2
    exit 2
fi

dataset=${datasets[$index]}
job_id=${SLURM_JOB_ID:-manual}
results_dir=${RSX_SLURM_RESULTS_DIR:-benchmarks/results/slurm}
mkdir -p "$results_dir"

# Default grids (override via env)
read -r -a priors <<< "${RSX_PRIORS:-0.001 0.005 0.01 0.02 0.05 0.1}"
read -r -a psexes <<< "${RSX_PSEX:-0.8 0.85 0.9 0.95}"

workdir=benchmarks/literature-workdir
popmap="$workdir/${dataset}/popmap.tsv"

# Prefer the precomputed d10-filtered table if it exists, otherwise fall back
# to markers_table.tsv (the analysis applies --min-depth 10 either way).
if [[ -f "$workdir/${dataset}/markers_d10.tsv" ]]; then
    markers="$workdir/${dataset}/markers_d10.tsv"
elif [[ -f "$workdir/${dataset}/markers_table.tsv" ]]; then
    markers="$workdir/${dataset}/markers_table.tsv"
else
    echo "Missing markers table for $dataset (neither markers_d10.tsv nor markers_table.tsv). Run literature_prepare_dataset.sbatch first." >&2
    exit 3
fi

if [[ ! -f "$popmap" ]]; then
    echo "Missing popmap for $dataset. Run literature_prepare_dataset.sbatch first." >&2
    exit 3
fi

for pi in "${priors[@]}"; do
  for p_sex in "${psexes[@]}"; do
    out_tsv="$results_dir/triage_${dataset}_pi${pi}_psex${p_sex}_${job_id}.tsv"
    echo "Running rsx triage for $dataset pi=$pi p_sex=$p_sex"

    "$pixi" run -e benchmark rsx triage \
      -t "$markers" \
      -p "$popmap" \
      -o "$out_tsv" \
      --min-depth 10 \
      --posterior-threshold 0.9 \
      --bayes-factor-threshold 10 \
      --prior-probability "$pi" \
      --linked-probability "$p_sex" \
      --groups M,F
  done
done

echo "Prior sensitivity runs for $dataset complete."
