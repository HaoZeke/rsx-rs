#!/usr/bin/env bash
#SBATCH --job-name=rsx-lit-prepare
#SBATCH --output=rsx-lit-prepare-%A_%a.out
#SBATCH --error=rsx-lit-prepare-%A_%a.err
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
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
    echo "SLURM_ARRAY_TASK_ID=$index is outside ${#datasets[@]} configured datasets" >&2
    exit 2
fi

dataset=${datasets[$index]}
job_id=${SLURM_JOB_ID:-manual}

if [[ ! -x target/release/rsx ]]; then
    echo "target/release/rsx not found; submit benchmarks/slurm/rsx_release_build.sbatch before this job" >&2
    exit 2
fi

"$pixi" run -e benchmark python benchmarks/run_literature_benchmarks.py \
    --manifest benchmarks/literature_datasets.tsv \
    --dataset "$dataset" \
    --workdir benchmarks/literature-workdir \
    --results "benchmarks/results/slurm/literature_benchmark_results_${dataset}_${job_id}.csv" \
    --summary "benchmarks/results/slurm/literature_benchmark_summary_${dataset}_${job_id}.tsv" \
    --force
