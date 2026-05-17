#!/usr/bin/env bash
#SBATCH --job-name=rsx-lit-radsex
#SBATCH --output=rsx-lit-radsex-%A_%a.out
#SBATCH --error=rsx-lit-radsex-%A_%a.err
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

: "${RADSEX_BIN:?Set RADSEX_BIN to the C++ RADSex executable path}"

read -r -a datasets <<< "${RSX_DATASETS:-notothenia_rossii danio_albolineatus tinca_tinca plecoglossus_altivelis}"
index=${SLURM_ARRAY_TASK_ID:-0}
if (( index < 0 || index >= ${#datasets[@]} )); then
    echo "SLURM_ARRAY_TASK_ID=$index is outside ${#datasets[@]} configured datasets" >&2
    exit 2
fi

dataset=${datasets[$index]}
job_id=${SLURM_JOB_ID:-manual}
results_dir=${RSX_SLURM_RESULTS_DIR:-benchmarks/results/slurm}
mkdir -p "$results_dir"
result="$results_dir/literature_speed_comparison_${dataset}_${job_id}.csv"

if [[ ! -x target/release/rsx ]]; then
    echo "target/release/rsx not found; submit benchmarks/slurm/rsx_release_build.sbatch before this job" >&2
    exit 2
fi

"$pixi" run -e benchmark python benchmarks/compare_literature_radsex.py \
    --workdir benchmarks/literature-workdir \
    --results "$result" \
    --rsx target/release/rsx \
    --radsex "$RADSEX_BIN" \
    --dataset "$dataset" \
    --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --force

echo "Wrote $result"
