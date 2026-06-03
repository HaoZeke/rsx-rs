# Minimal demo for rsx / pyrsx

This directory contains (or generates) a tiny synthetic RAD-seq-like dataset so you can try the full pipeline in < 30s without downloading real data.

## CLI demo

```bash
# 1. Generate tiny data (or use the committed tiny/ if present)
python ../../benchmarks/generate_data.py --outdir data --scale small --nind 8 --nmarkers 5000

# 2. Build marker table
rsx process -i data/reads -o markers.tsv -T 4 -d 2

# 3. Basic stats + sex signal
rsx freq -t markers.tsv -o freq.tsv -d 2
echo -e "ind1\tM\nind2\tM\nind3\tM\nind4\tM\nind5\tF\nind6\tF\nind7\tF\nind8\tF" > popmap.tsv
rsx distrib -t markers.tsv -p popmap.tsv -o distrib.tsv -G M,F -d 2
rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv -G M,F -d 2 --bayes

# 4. (optional) triage
rsx triage -t markers.tsv -p popmap.tsv -o triage.tsv -G M,F
```

## Python high-level demo

```python
import pyrsx
from pathlib import Path

pyrsx.process("data/reads", "markers.tsv", threads=4, min_depth=2)

tbl = pyrsx.MarkerTable.from_path("markers.tsv")
print("markers:", len(tbl))

# distrib / signif with bayes
pyrsx.distrib("markers.tsv", "popmap.tsv", "distrib.tsv", groups=["M","F"])
pyrsx.signif("markers.tsv", "popmap.tsv", "signif.tsv", groups=["M","F"], bayes=True)

# PCA (streaming Tucker)
pyrsx.pca("markers.tsv", "pca_out/", n_components=3)
```

See the main quickstart and https://rsx.rgoswami.me for the full story and memory guarantees.
