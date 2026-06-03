#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

echo "== rsx minimal demo =="
mkdir -p data

if command -v rsx >/dev/null; then
    echo "Using rsx from PATH"
else
    echo "rsx not on PATH; build first (cargo build --release -p rsx-cli) or install"
    exit 1
fi

python ../../benchmarks/generate_data.py --outdir data --scale small --nind 8 --nmarkers 5000 || {
    echo "generate_data.py not found or failed; using pre-existing tiny data if any"
}

rsx process -i data/reads -o markers.tsv -T 4 -d 2
echo -e "ind1\tM\nind2\tM\nind3\tM\nind4\tM\nind5\tF\nind6\tF\nind7\tF\nind8\tF" > popmap.tsv
rsx distrib -t markers.tsv -p popmap.tsv -o distrib.tsv -G M,F -d 2
rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv -G M,F -d 2 --bayes

echo "CLI done. markers=$(wc -l < markers.tsv) lines"

if python -c 'import pyrsx' 2>/dev/null; then
    python - << 'PY'
import pyrsx
pyrsx.pca("markers.tsv", "pca_demo/", n_components=2)
print("Python PCA demo wrote pca_demo/")
PY
fi

echo "Minimal demo complete. See README.md here and the main docs."
