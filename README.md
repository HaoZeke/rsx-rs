<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/branding/rsx-logo-dark.svg">
    <img src="assets/branding/rsx-logo.svg" alt="rsx" width="460">
  </picture>
</p>

**High-performance streaming toolkit for RAD-seq sex determination.**

A Rust framework for RAD-seq marker analysis and sex determination: bounded-memory streaming kernels, a Bayesian marker-evidence layer, and Python and C bindings. Builds on and stays command-compatible with [RADSex](https://github.com/RomainFeron/RADSex), so prior results remain directly comparable.

[![CI](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/HaoZeke/rsx-rs/graph/badge.svg)](https://app.codecov.io/gh/HaoZeke/rsx-rs)
[![Documentation](https://img.shields.io/badge/docs-rsx.rgoswami.me-blue)](https://rsx.rgoswami.me)
[![Crates.io](https://img.shields.io/crates/v/rsx-cli?label=crates.io)](https://crates.io/crates/rsx-cli)
[![PyPI](https://img.shields.io/pypi/v/pyrsx)](https://pypi.org/project/pyrsx/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Install

### From GitHub releases (recommended for end users)
Pre-built binaries for Linux (x86_64/aarch64), macOS (x86_64/arm64), Windows (without `map`):

```bash
# See https://github.com/HaoZeke/rsx-rs/releases for the latest
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/HaoZeke/rsx-rs/releases/download/v0.2.6/rsx-cli-installer.sh | sh
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

See the [Python README](rsx-python/README.md) for the `MarkerTable`,
`TriageResult`, narwhals, and plotting APIs.

## 30-second quickstart

```bash
# CLI
rsx process -i reads/ -o markers.tsv -T 8 -d 5
rsx distrib -t markers.tsv -p popmap.tsv -o distrib.tsv -G M,F
rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv -G M,F --bayes
rsx map -t markers.tsv -p popmap.tsv -g genome.fa -o aligned.tsv -G M,F

# Optional NVIDIA CUDA path for chi-square significance testing
cargo build -p rsx-cli --release --features cuda
target/release/rsx signif -t markers.tsv -p popmap.tsv -o signif.tsv \
  -G M,F --test chisq --backend cuda

# Five measured repetitions across the GPU crossover range
cargo run -p rsxcore --release --example benchmark_compute_backends \
  --features cuda -- 1000,10000,100000,1000000,10000000 5

# Python
import pyrsx
pyrsx.process("reads/", "markers.tsv", threads=8, min_depth=5)
pyrsx.signif("markers.tsv", "popmap.tsv", "signif.tsv", test="fisher", correction="fdr", bayes=True)
tbl = pyrsx.MarkerTable.from_path("markers.tsv")
tbl.summary()
```

Full pipeline, memory guarantees, and all 10 commands (including new `merge`, `pca`, `triage`) are documented at https://rsx.rgoswami.me .

## Versioned run profiles

Every analysis argument fits in a strict, schema-versioned Tom's Obvious
Minimal Language (TOML) profile. Explicit command-line values override the
profile. The resolved profile records the values that the calculation
receives.

```bash
rsx --profile examples/profiles/triage.toml \
  --write-hydrated-profile triage.hydrated.toml \
  --reproducibility-archive triage-reproducibility.zip
```

The hydrated TOML includes compatibility defaults such as all
Bayesian prior shapes and the depth streaming policy. For example, `depth`
accepts `--streaming-mode auto|memory|streaming` and
`--streaming-threshold-bytes`; both values are preserved in a hydrated
profile.

rsx creates the reproducibility ZIP before input validation or analysis.
A failed invocation therefore retains its resolved configuration whenever
profile resolution succeeds. The ZIP contains the input and hydrated TOML,
executable, lockfile, build and run manifests, JSON Schema, CycloneDX
software bill of materials, citation and license files, and checksums. It
omits analysis results and input datasets.

## Features

- All original RADSex commands + `merge` (external sort for 75M+ markers, ~500 MB RAM), `pca` (streaming sample-space / Tucker mode-2 factors), `triage` (Bayes + strict candidate ranking).
- Bounded-memory streaming for the default analysis paths (Bonferroni `signif`, `distrib`, `freq`, `triage`, streaming `pca` / `merge`). Exceptions: FDR `signif` stores O(n_markers) p-values for BH correction then re-streams; `depth` can use a configurable automatic threshold, forced memory mode, or forced streaming mode. See the [quickstart memory table](https://rsx.rgoswami.me/tutorials/quickstart.html).
- More than twofold faster than C++ RADSex on literature panels, with matching output when groups are specified.
- Python (`pyrsx`), R (`rsxr`), and C (`rsx.h` / cargo-c) bindings over the shared core.
- Optional: NVIDIA CUDA chi-square batches for `signif`, parquet I/O, MPI,
  and minimap2 mapping (controlled by feature flags on Windows). CUDA selection is explicit;
  unsupported tests and builds without the `cuda` feature return an error. The
  CUDA path transfers count pairs directly, retains the compiled kernel, and
  reuses the largest page-locked result buffer within a process. Metrics keep
  first-batch setup separate from transfers and kernel execution.
- Reproducible: pixi environments, an Airspeed Velocity literature benchmark suite, and SymPy derivations and a Lean proof for the math.

## Documentation

- Full site: https://rsx.rgoswami.me (tutorials, command reference, architecture, HPC design, R + Python integration).
- Preprint: Goswami R, Goswami R. rsx: A high-performance streaming toolkit for RAD-seq sex determination. arXiv:2606.06434 (2026). https://arxiv.org/abs/2606.06434 (submitted to BMC Bioinformatics).
- Reproducibility materials: the companion `rsx_bmc_repro` package (snakemake-orchestrated, MCA/Zenodo archive shape matching the rest of the collection) + the org files under `repro/` in this repo.
- Software archive: https://doi.org/10.5281/zenodo.20531538

## Ecosystem

rsx is the RAD-seq sex-marker engine of a salmonid sex-determination toolkit.
Two R packages consume its outputs and cover the downstream analyses:

- [ChromSex](https://github.com/RuhiRG/ChromSex) -- sex-balanced Hudson / Weir & Cockerham FST, sex-specific heterozygosity, and genomic content around sex loci.
- [SexPCR](https://github.com/RuhiRG/SexPCR) -- sdY PCR primer screen and assay audit against dissection sex.

The brown trout LG28 manuscript applies all three on one Icelandic dataset.

## Citation

Please cite the preprint (or published version when available).

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

For the software itself (v0.2.6), prefer the GitHub Cite button, the `CITATION.cff`, or the Zenodo DOI entry (generated from the deposit page). The reproducibility archive has its own Zenodo-generated BibTeX.

## Contributing

See `CONTRIBUTING.md`. We use pixi for dev envs, conventional commits, and the usual Rust `cargo fmt && cargo clippy -D warnings && cargo test`.

## License

GPL-3.0-or-later. See `LICENSE`.

The C++/Python/RADSex heritage is similarly licensed; see original sources.
