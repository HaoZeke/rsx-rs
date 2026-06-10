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
curl -sSfL https://github.com/HaoZeke/rsx-rs/releases/download/v0.2.3/rsx-installer.sh | sh
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
- Preprint: Goswami R, Goswami R. rsx: A high-performance streaming toolkit for RAD-seq sex determination. arXiv:2606.06434 (2026). https://arxiv.org/abs/2606.06434 (submitted to BMC Bioinformatics).
- Reproducibility materials: the companion `rsx_bmc_repro` package (snakemake-orchestrated, MCA/Zenodo archive shape matching the rest of the collection) + the org files under `repro/` in this repo.
- Software archive: https://doi.org/10.5281/zenodo.20531538

## Citation

Please cite the preprint (or published version when available) and the RADSex reference where relevant.

> Goswami R, Goswami R. /rsx: A high-performance streaming toolkit for RAD-seq
> sex determination./ arXiv:2606.06434 (2026). https://arxiv.org/abs/2606.06434

See `CITATION.cff` (root) for the machine-readable entry (includes the arXiv preprint, RADSex, and Zenodo software DOI). 
GitHub's "Cite this repository" button also generates BibTeX/APA from it.

RADSex reference: Feron et al., Mol Ecol Resour 2021. https://doi.org/10.1111/1755-0998.13360

For the benchmark data / figures, also cite the deposited reproducibility archive (Zenodo): https://doi.org/10.5281/zenodo.20531539 .

### BibTeX / BibLaTeX

```bibtex
@article{Goswami2026rsx,
  title         = {rsx: A high-performance streaming toolkit for RAD-seq sex determination},
  author        = {Goswami, Rohit and Goswami, Ruhila},
  year          = {2026},
  eprint        = {2606.06434},
  archivePrefix = {arXiv},
  primaryClass  = {q-bio.GN},
  url           = {https://arxiv.org/abs/2606.06434},
  doi           = {10.48550/arXiv.2606.06434},
  note          = {Preprint, submitted to BMC Bioinformatics}
}
```

For the software itself (v0.2.3), prefer the GitHub Cite button, the `CITATION.cff`, or the Zenodo DOI entry (generated from the deposit page). The reproducibility archive has its own Zenodo-generated BibTeX.

## Contributing

See `CONTRIBUTING.md`. We use pixi for dev envs, conventional commits, and the usual Rust `cargo fmt && cargo clippy -D warnings && cargo test`.

## License

GPL-3.0-or-later. See `LICENSE`.

The C++/Python/RADSex heritage is similarly licensed; see original sources.
