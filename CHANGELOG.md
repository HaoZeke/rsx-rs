# Changelog

All notable changes to rsx-rs are documented here.

## [0.2.9] - 2026-08-02

### Added
- **CUDA / signif**: optional `cuda` feature and explicit
  `signif --backend cuda` path for batched `f64` Yates chi-square p-values;
  record setup, transfer, kernel, and return timings without silent fallback.

### Fixed
- **depth**: exact and streaming modes report the mathematical median for even
  sample counts by averaging the two middle depths; half-integer medians are
  preserved in TSV output.
- **triage**: preserve the strongest marker as an explicitly non-call
  `exploratory` row when every marker falls below the configured strict,
  posterior, and Bayes-factor thresholds.

### Removed


- The Sollya-fitted degree-40 erfc polynomial, its coefficient table, the
  Horner helper, the `hexf` dependency, and `scripts/sollya/`. The bound the
  coefficients shipped with, 8.22e-17, describes the mathematical polynomial;
  evaluated in binary64 Horner form over the fitted [0, 6] the routine left that
  bound at t ~ 0.36 and returned negative values from t ~ 4, which a probability
  cannot be. Nothing called it -- `chi_squared_p` goes through libm and the CUDA
  kernels call the device `erfc` -- so no result changes.

### Changed


- `tests/test_precision.rs` asserts the shipped path tracks `erfc` to 1e-15 and
  stays inside [0, 1] across the same range the polynomial was fitted on.

## [0.2.6] - 2026-07-06

Correctness and robustness patch (process table metadata and error handling,
`signif` thresholds/FDR streaming, C FFI UTF-8 validation, write-error
propagation, `distrib` flag parity) plus full-suite nextest in CI. Software
Zenodo concept DOI unchanged (`10.5281/zenodo.20531538`); deposit a new
**version** of that record for 0.2.6. Whether the reproducibility archive
(`10.5281/zenodo.20531539`) needs a new version depends on which workflows
you re-run (see Fixed below — process/`min_depth` and FDR paths can change
outputs).

### Fixed
- **process**: `#Number of markers` header matches rows retained after
  `min_depth`; fail the run if any individual input file errors; deterministic
  row order by packed sequence key.
- **signif**: compare p-values in `f64` (no `f32` cast); FDR stores p-values
  only then re-streams (two table passes); reject FASTA+FDR *before* truncating
  the output path; propagate write errors.
- **C API**: reject non-UTF-8 C strings with `RSX_INVALID_PARAMETER` (never map
  to empty path/group).
- **stats**: Fisher documented and tested as two-sided (density method);
  `g_test`/`fisher_exact` return `p=1` on invalid tables; `group_bias` safe on
  empty groups.
- **triage / subset / map**: propagate output write failures instead of
  discarding `io::Result`.

### Added
- **distrib** CLI: `--test`, `--correction`, `--bayes` (parity with Python/R);
  cell-grid FDR and optional Bayes columns; document that Bonferroni still uses
  `n_markers` (RADSex-compatible) while FDR uses the cell grid.
- **pca**: also emit `sample_scores.tsv` (same matrix as legacy `loadings.tsv`).
- **CI / dev**: `cargo-nextest` wired in pixi and GitHub Actions for the full
  `rsxcore` suite (lib + pipeline + precision); `.config/nextest.toml`.
- Docs site: `sphinxcontrib-bibtex` + `docs/source/references.bib`; README
  memory claims and architecture mmap narrative corrected.

### Changed
- Integration tests use unique `tempfile` directories (nextest-safe).
- k-mer dedup wording: min-hash LSH, not “majority”.

## [0.2.5] - 2026-07-01

Corrective lockstep cut superseding the broken `0.2.4` publish. crates.io
cannot replace a version once uploaded, so `rsxcore`/`rsx-cli` `0.2.4` are
yanked and this release is the CI-clean source of truth.

### Changed
- **Median**: exact upper-median via `select_nth_unstable_by` (no full sort of
  the depth vector); order-stat identity checked with Lean and an independent
  SymPy/quickselect validator (`proofs/lean/MedianSelect`,
  `scripts/sympy/median_select_proof.py`).
- **Fisher / signif**: online log-sum-exp accumulation (single pass, no second
  materialised log-prob buffer).
- **TSV parse**: `fast_parse_u16` for depth cells on the hot path.
- **Depth**: min/max without an extra full sort before the median selection.
- **PCA**: streaming Gram contraction (`accumulate_gram_streaming`) with a
  sparse rank-1 microbench path; dense panel path kept for typical RAD panels.
- **Bitset**: unrolled popcount on masked words for triage/signif masks.
- **Frequency (Arrow)**: columnar tile histogram path for marker streams.
- **Packaging**: `rsx-cli` depends on `rsxcore` `0.2.5` in lockstep (the
  yanked `rsx-cli` `0.2.4` on crates.io still required `rsxcore` `^0.2.3`
  while `rsxcore` `0.2.4` had already been published from a failing commit).

### Fixed
- Sphinx/Shibuya docs: code-fence whitespace (unwrap `data-line` / force
  inline tokens), single frame on `pre` blocks (override theme `display:grid`),
  and readable token colours on dark backgrounds (custom.css after doxyrest;
  inherit doxyrest keyword styles).
- Docs linkcheck (lychee) and org→RST conversion issues for the published site.
- `cargo fmt` and `clippy -D warnings` cleanups that landed after the
  mistaken `v0.2.4` tag (`checked_div` in depth averages, infallible `Ok`
  return in frequency, PCA doc/`GramStreamingAcc` type alias, Arrow source
  formatting).

### Added
- rOpenSci-oriented **rsxr** scaffolding and monorepo-aware pkgcheck CI (package
  remains at 0.1.0 pending software-review submission).
- Language-bindings matrix and readcon-core-aligned docs flourish on the site.

## [0.2.4] - 2026-07-01

YANKED on crates.io and superseded by `0.2.5`. The `v0.2.4` tag pointed at
`db230c6` (`chore(release): 0.2.4`) before the `rustfmt`/`clippy` fixes on
`main`, so CI never passed on that commit; `rsxcore` `0.2.4` was published
from it while `rsx-cli` `0.2.4` still depended on `rsxcore` `^0.2.3`.

## [0.2.3] - 2026-06-04

### Added
- Accompanying preprint: Goswami R, Goswami R. rsx: A high-performance
  streaming toolkit for RAD-seq sex determination. arXiv:2606.06434 (2026).
  https://arxiv.org/abs/2606.06434 (submitted to BMC Bioinformatics).
  Software Zenodo: https://doi.org/10.5281/zenodo.20531538 .
  Reproducibility archive: https://doi.org/10.5281/zenodo.20531539 .

### Fixed
- Prior-sensitivity heatmap now picks the cell annotation colour by value
  (white on the dark high-count cells, near-black on light cells) so the
  counts stay legible across the full fill range.

## [0.2.2] - 2026-06-04

### Fixed
- Retained the `depth` command rows in the literature speedup summary and pairs
  tables, so the reproducibility workflow regenerates all 56 paired command/
  dataset/depth timings reported in the paper instead of dropping the `depth`
  rows.

## [0.2.1] - 2026-06-03

### Fixed
- Hardened TSV result handling so Python bindings preserve core command headers
  and load core TSV outputs through Arrow-compatible readers.
- Guarded Python extension rebuilds against stale artifacts when switching Pixi
  environments.
- Aligned the Python package metadata with the `0.2.1` release version.
- Rejected malformed `--groups` CLI values instead of silently falling back to
  popmap group auto-detection.
- Made benchmark manifests tolerant of comments and added ENA FASTQ retry logic
  for literature benchmark downloads.

### Added
- Reproducibility package validation for the tracked benchmark payloads,
  including the Python build helper needed by archive users.
- Focused regression coverage for source parity fixtures and Python benchmark
  package layout.

## [0.2.0] - 2026-06-03

### Added
- **`rsx pca` command**: streaming PCA/Tucker mode-2 decomposition of the
  depth matrix. Computes Gram matrix X^T X streaming (O(n_ind^2) memory),
  eigendecomposes for exact Tucker mode-2 factors. Use cases: sex signal
  detection, sample QC, population structure analysis.
  SymPy proof: `scripts/sympy/tucker_covariance_proof.py`.
- **Sparse external sort for depth**: only non-zero depths sorted (zeros
  tracked by count), 3.3x I/O reduction for 70% sparse RAD-seq data.
  SymPy proof: `scripts/sympy/sparse_median_proof.py`.
- **Bounded-memory streaming for ALL commands.** No command accumulates
  O(n_markers) data in memory regardless of input size.
  - signif/subset: two-pass streaming with Bonferroni count in pass 1
  - map: two-pass (count then align+write)
  - depth: exact median via external sort for files > 2GB, auto-detected
  - merge: external sort with lz4 temp files + k-way merge
- External sort-merge for `merge` command: bounded-memory (~500MB) merge
  of 75M+ unique sequences using chunked sort + lz4 temp files + k-way merge.
  Fixes OOM on large RAD-seq datasets (25GB+ input).
- `--buffer-size` flag for merge command to tune memory/temp-file tradeoff.
- `--output-parquet` flag for merge: Parquet output with ZSTD compression
  (feature-gated behind `parquet-io`).
- 2-bit DNA packing: 100bp sequences stored as 26 bytes (4x compression).
  Used in both `process` and `merge` commands.
- DashMap concurrent merge in `process` command (for >= 8 individuals).
- Optional MPI support for `process` command (`--features mpi`).
- Feature-gated minimap2: `map` feature (default on) allows Windows builds
  without minimap2 for all other commands.
- Workspace package version inheritance in root Cargo.toml for consistent
  versioning across crates.

### Changed
- Merge command input files are now positional arguments (glob-friendly):
  `rsx merge -o out.tsv file1.tsv file2.tsv` instead of `-i file1 file2`.
- Merge stores depths as `u16` instead of `String` (memory reduction).
- signif/subset are ~30% slower on large data due to two-pass overhead,
  but use O(n_individuals) memory instead of O(n_markers).
- Python package version aligned to 0.2.0 (was 0.1.0); Sphinx release and
  all Cargo crates now consistently 0.2.0 via workspace.package.

## [0.1.0] - 2026-04-06

### Added
- All 7 RADSex commands: process, distrib, signif, depth, freq, map, subset.
- Byte-identical TSV output to C++ RADSex (when groups specified explicitly).
- 2.6x faster than C++ RADSex across all benchmarks.
- Bitset + popcount group counting (no HashMap in hot path).
- mmap zero-copy I/O for markers table parsing.
- SymPy-derived erfc identity for chi-squared CDF (eliminates gamma function).
- Sollya degree-40 minimax polynomial for GPU erfc (via `hexf!` macro).
- Three-tier parsing: fast (min_depth=1, skip id+seq), medium, full.
- minimap2 alignment (replaces BWA-MEM).
- Parallel file processing via rayon.
- C API via cbindgen for R/Python/C++ integration.
- Cross-platform releases: Linux (x86_64, aarch64), macOS (x86_64, arm64).
- ASV benchmarks with CI regression detection.
- 49 tests (31 unit + 9 integration + 9 precision).
- Orgmode documentation with Sphinx export.
- SymPy derivation scripts and Sollya minimax scripts.

### Performance (vs C++ RADSex v1.2.0)
- Overall: 2.6x faster
- distrib: 4.0x
- depth: 4.2x
- process: 2.8x (rayon parallel)
- map: 2.1x (minimap2)
