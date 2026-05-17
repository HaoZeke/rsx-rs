#!/usr/bin/env bash
#SBATCH --job-name=rsx-release-build
#SBATCH --output=rsx-release-build-%j.out
#SBATCH --error=rsx-release-build-%j.err
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
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

"$pixi" run -e benchmark build-rsx
