#!/usr/bin/env bash
# literature_biology_low_depth.sbatch
# Run the full literature analysis pipeline at low depths (3,5,8,10) on real data.
# Produces data for low-depth power curves.
#
#SBATCH --job-name=rsx-bio-lowdepth
#SBATCH --output=rsx-bio-lowdepth-%A_%a.out
#SBATCH --error=rsx-bio-lowdepth-%A_%a.err
#SBATCH --partition=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
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
    echo "pixi not found" >&2
    exit 127
fi

read -r -a datasets <<< "${RSX_DATASETS:-danio_albolineatus notothenia_rossii plecoglossus_altivelis tinca_tinca}"
index=${SLURM_ARRAY_TASK_ID:-0}
dataset=${datasets[$index]}

read -r -a depths <<< "${RSX_DEPTHS:-3 5 8 10}"

workdir=benchmarks/literature-workdir

"$pixi" run -e paper build-python

for d in "${depths[@]}"; do
  echo "=== $dataset min_depth=$d ==="

  # Re-use the existing mode analysis but with custom min_depth
  "$pixi" run -e paper python benchmarks/analyze_literature_modes.py \
    --workdir "$workdir" \
    --results "benchmarks/results/lowdepth_${dataset}_d${d}.csv" \
    --dataset "$dataset" \
    --min-depth "$d"
done

echo "Low-depth runs for $dataset complete."
