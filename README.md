<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/branding/rsx-logo-dark.svg">
    <img src="assets/branding/rsx-logo.svg" alt="rsx" width="460">
  </picture>
</p>

**High-performance streaming toolkit for RAD-seq sex determination.**

A Rust framework for RAD-seq marker analysis and sex determination: bounded-memory streaming kernels, a Bayesian marker-evidence layer, and Python and C bindings. Builds on and stays command-compatible with [RADSex](https://github.com/RomainFeron/RADSex), so prior results remain directly comparable.

[![CI](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-rsx.rgoswami.me-blue)](https://rsx.rgoswami.me)
[![Crates.io](https://img.shields.io/crates/v/rsx-cli?label=crates.io)](https://crates.io/crates/rsx-cli)
[![PyPI](https://img.shields.io/pypi/v/pyrsx)](https://pypi.org/project/pyrsx/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Install

### From GitHub releases (recommended for end users)
Pre-built binaries for Linux (x86_64/aarch64), macOS (x86_64/arm64), Windows (without `map`):

```bash
# See https://github.com/HaoZeke/rsx-rs/releases for the latest
curl -sSfL https://github.com/HaoZeke/rsx-rs/releases/download/v0.2.1/rsx-installer.sh | sh
```

### From source

```bash
git clone https://github.com/HaoZeke/rsx-rs.git
cd rsx-rs
cargo build --release
# binary at target/release/rsx
```

### Via pixi (dev / reproducible)

```bash
pixi run build
# or pixi run -e dev build-portable
```

### Python bindings

```bash
pip install pyrsx
```

See the [Python README](rsx-python/README.md) for the high-level `MarkerTable`, `TriageResult`, narwhals, and plotting APIs.

## 30-second quickstart

```bash
# CLI
rsx process -i reads/ -o markers.tsv -T 8 -d 5
rsx distrib -t markers.tsv -p popmap.tsv -o distrib.tsv -G M,F
rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv -G M,F --bayes
rsx map -t markers.tsv -p popmap.tsv -g genome.fa -o aligned.tsv -G M,F

# Python (high-level)
import pyrsx
pyrsx.process("reads/", "markers.tsv", threads=8, min_depth=5)
pyrsx.signif("markers.tsv", "popmap.tsv", "signif.tsv", test="fisher", correction="fdr", bayes=True)
tbl = pyrsx.MarkerTable.from_path("markers.tsv")
...
```

Full pipeline, memory guarantees, and all 10 commands (including new `merge`, `pca`, `triage`) are documented at https://rsx.rgoswami.me .

## Features

- All original RADSex commands + `merge` (external sort for 75M+ markers, ~500 MB RAM), `pca` (streaming Tucker), `triage` (Bayes + strict candidate ranking).
- Bounded-memory streaming for every command — no O(n_markers) accumulation.
- 2-6x+ faster than C++ RADSex on literature panels (byte-identical output when groups specified).
- Python bindings (low-level + ergonomic `MarkerTable` / Arrow / narwhals), C API via cbindgen.
- Optional: parquet I/O, MPI, minimap2 mapping (feature-gated for Windows).
- Reproducible: pixi environments, ASV + literature benchmark harness, SymPy/Sollya proofs for the math.

## Documentation

- Full site: https://rsx.rgoswami.me (tutorials, command reference, architecture, HPC design, R + Python integration).
- Paper (BMC Bioinformatics, in submission): see the companion manuscript repository or the forthcoming published version.
- Reproducibility materials: the companion `rsx_bmc_repro` package (snakemake-orchestrated, MCA/Zenodo archive shape matching the rest of the collection) + the org files under `repro/` in this repo.

## Citation

Please cite the software article (when published) and the RADSex reference it extends.

See `CITATION.cff` (root) for the machine-readable entry.

For the benchmark data / figures used in the paper, also cite the deposited reproducibility archive (Zenodo / Materials Cloud Archive entry, to be minted from the `rsx_bmc_repro` package after the heavy builder runs).

RADSex reference: Feron et al., Mol Ecol Resour 2021. https://doi.org/10.1111/1755-0998.13360

## Contributing

See `CONTRIBUTING.md`. We use pixi for dev envs, conventional commits, and the usual Rust `cargo fmt && cargo clippy -D warnings && cargo test`.

## License

GPL-3.0-or-later. See `LICENSE`.

The C++/Python/RADSex heritage is similarly licensed; see original sources.
