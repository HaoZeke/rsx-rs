# pyrsx

Python bindings for [rsx](https://github.com/HaoZeke/rsx-rs): a high-performance streaming toolkit for RAD-seq sex determination.

## Installation

```bash
pip install pyrsx
```

## Usage

```python
import pyrsx

# Process FASTQ files into marker depth table
pyrsx.process("reads/", "markers.tsv", threads=4, min_depth=5)

# Compute distribution with Fisher's exact test + FDR
pyrsx.distrib("markers.tsv", "popmap.tsv", "distrib.tsv",
              test="fisher", correction="fdr")

# Extract significant markers with Bayesian output
pyrsx.signif("markers.tsv", "popmap.tsv", "signif.tsv",
             test="fisher", correction="fdr", bayes=True)

# Streaming PCA
pyrsx.pca("markers.tsv", "pca_results/", n_components=10)

# Merge tables (bounded memory, handles 75M+ sequences)
pyrsx.merge(["table1.tsv", "table2.tsv"], "merged.tsv")
```

## Features

- All rsx commands accessible from Python
- 3.14x geometric-mean speedup on the tracked Slurm literature comparison panel
- Bounded-memory streaming for arbitrarily large datasets
- Multiple statistical tests: chi-squared, Fisher's exact, G-test
- Multiple corrections: Bonferroni, Benjamini-Hochberg FDR
- Bayesian sex-linkage classification (Bayes Factor + posterior)
- Streaming PCA via Tucker mode-2 decomposition
- K-mer based marker deduplication
