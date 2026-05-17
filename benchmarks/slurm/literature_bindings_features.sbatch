#!/usr/bin/env bash
#SBATCH --job-name=rsx-lit-bindings
#SBATCH --output=rsx-lit-bindings-%j.out
#SBATCH --error=rsx-lit-bindings-%j.err
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
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

dataset=${RSX_DATASET:-notothenia_rossii}
min_depth=${RSX_MIN_DEPTH:-10}
n_components=${RSX_N_COMPONENTS:-4}

"$pixi" run -e paper build-python
"$pixi" run -e paper python benchmarks/run_literature_bindings.py \
    --workdir benchmarks/literature-workdir \
    --results "benchmarks/results/literature_binding_results_${dataset}_${SLURM_JOB_ID:-manual}.csv" \
    --output-dir benchmarks/literature-workdir/bindings \
    --dataset "$dataset" \
    --min-depth "$min_depth" \
    --n-components "$n_components"
