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

## High-level API & backend agnosticism (recommended)

The low-level functions above are thin wrappers. For most users the
`MarkerTable` + result objects (in `pyrsx.api`) are the idiomatic entry
point:

```python
import pyrsx as rsx

table = rsx.MarkerTable.from_path("markers.tsv")   # or from_dataframe(...)
result = table.triage(popmap="popmap.tsv", min_depth=10)

# Everything is a narwhals DataFrame under the hood → backend agnostic
print(result.df)                    # stays in whatever backend you prefer
df = result.to_polars()             # or .to_pandas(), to_dataframe(backend=...)
```

**How outputs are read (no forced pandas fallback):**
Internal TSVs produced by rsx core commands are read with pyarrow.csv
(handling the leading `#Number of markers` comment via `skip_rows=1`)
and then wrapped with `to_narwhals(...)`. The exposed objects are
always narwhals DataFrames (concrete backend = pyarrow by default for
efficiency). You only pull in pandas/polars if *you* ask for that
backend later. This is the standard narwhals approach used throughout
the high-level API (see `_adapters.py`, `_read_core_tsv`, and the
detailed docs in the Rust extension).

See the docstrings of `MarkerTable`, the various `*Result` classes, and
`_read_core_tsv` for the full rationale.
