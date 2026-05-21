#!/usr/bin/env bash
#SBATCH --job-name=rsx-lit-bindings
#SBATCH --output=rsx-lit-bindings-%A_%a.out
#SBATCH --error=rsx-lit-bindings-%A_%a.err
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
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

if [[ -n "${RSX_DATASET:-}" ]]; then
    datasets=("$RSX_DATASET")
else
    read -r -a datasets <<< "${RSX_DATASETS:-danio_albolineatus notothenia_rossii plecoglossus_altivelis tinca_tinca}"
fi
index=${SLURM_ARRAY_TASK_ID:-0}
if (( index < 0 || index >= ${#datasets[@]} )); then
    echo "SLURM_ARRAY_TASK_ID=$index is outside ${#datasets[@]} configured datasets" >&2
    exit 2
fi

dataset=${datasets[$index]}
job_id=${SLURM_JOB_ID:-manual}
results_dir=${RSX_SLURM_RESULTS_DIR:-benchmarks/results/slurm}
min_depth=${RSX_MIN_DEPTH:-10}
n_components=${RSX_N_COMPONENTS:-4}
mkdir -p "$results_dir"

"$pixi" run -e paper build-python
"$pixi" run -e paper python benchmarks/run_literature_bindings.py \
    --workdir benchmarks/literature-workdir \
    --results "$results_dir/literature_binding_results_${dataset}_${job_id}.csv" \
    --output-dir benchmarks/literature-workdir/bindings \
    --dataset "$dataset" \
    --min-depth "$min_depth" \
    --n-components "$n_components"
