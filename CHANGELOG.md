# Changelog

All notable changes to rsx-rs are documented here.

## [Unreleased]

### Added
- External sort-merge for `merge` command: bounded-memory (~500MB) merge
  of 75M+ unique sequences using chunked sort + lz4 temp files + k-way merge.
  Fixes OOM on large ChromSex datasets (25GB+ input).
- `--buffer-size` flag for merge command to tune memory/temp-file tradeoff.
- 2-bit DNA packing: 100bp sequences stored as 26 bytes (4x compression).
  Used in both `process` and `merge` commands.
- DashMap concurrent merge in `process` command (for >= 8 individuals).
  Rayon threads insert directly into shared map, eliminating sequential
  merge bottleneck.
- Optional MPI support for `process` command (`--features mpi`).
  Distributes FASTQ file processing across MPI ranks.

### Changed
- Merge command input files are now positional arguments (glob-friendly):
  `rsx merge -o out.tsv file1.tsv file2.tsv` instead of `-i file1 file2`.
- Merge stores depths as `u16` instead of `String` (memory reduction).

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
