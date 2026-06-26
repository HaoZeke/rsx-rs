# Introduction to rsxr

## Statement of need

RAD-seq sex-determination analyses need marker depth tables, group
contrasts, and significance or Bayesian triage. The
[rsx](https://github.com/HaoZeke/rsx-rs) engine implements those steps
as streaming, bounded-memory commands with a stable C API. **rsxr**
binds that C API directly from R (no subprocess), so pipelines can call
`process`, `freq`, `distrib`, `signif`, `triage`, `depth`, `merge`, and
`pca` without leaving R, and consume results as tibbles via
[`marker_table()`](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md).

Target users are R-first bioinformaticians who already work with popmaps
and marker TSVs (including outputs compatible with the classic RADSex
workflow) and want an in-process R entry point to the Rust
implementation.

## Installation

rsxr lives in the `rsx-r/` subdirectory of the rsx-rs monorepo and
requires a Rust toolchain (`cargo`, `rustc`; see <https://rustup.rs>)
plus the `rsxcore` sources as a sibling of the package directory
(default in a full clone).

``` r

# pak (recommended)
pak::pak("HaoZeke/rsx-rs/rsx-r")

# remotes
remotes::install_github("HaoZeke/rsx-rs", subdir = "rsx-r")
```

From a local clone of the monorepo:

``` bash
cd rsx-r
R CMD INSTALL .
# or: pixi run install
```

Set `RSX_CORE_DIR` if `rsxcore` is not at `../rsxcore` relative to the
package.

## Low-level API

Each function mirrors one rsx CLI command and returns the output path
invisibly:

``` r

library(rsxr)

rsx_version()

# After demultiplexed FASTQ files exist under reads/:
rsx_process("reads/", "markers.tsv", threads = 8L, min_depth = 5L)
rsx_freq("markers.tsv", "freq.tsv", min_depth = 5L)
rsx_signif(
  "markers.tsv", "popmap.tsv", "signif.tsv",
  test = "fisher", correction = "fdr", bayes = TRUE
)
rsx_triage("markers.tsv", "popmap.tsv", "triage.tsv", min_depth = 10L)
```

## High-level `marker_table` workflow

[`marker_table()`](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md)
holds a path (or writes a data frame to a temp TSV). Verbs run the
corresponding command into a temporary file and return a tibble:

``` r

mt <- marker_table("markers.tsv")

triaged <- triage(mt, popmap = "popmap.tsv", min_depth = 10L)
sig     <- signif_markers(mt, popmap = "popmap.tsv", test = "fisher")
dist    <- distrib(mt, popmap = "popmap.tsv", group1 = "M", group2 = "F")
freqs   <- frequencies(mt, min_depth = 5L)
depths  <- depth(mt, popmap = "popmap.tsv")
```

## System requirements and CRAN builds

- **Runtime build:** `cargo` and `rustc` (see `Config/rsxr/MSRV` in
  `DESCRIPTION`).
- **Compression libs** used by the Rust stack: zlib, bzip2, xz (usually
  present on CRAN builders; install `-dev` packages on Linux if linking
  fails).
- **Release tarballs** may ship vendored Rust crates under
  `src/rust/vendor.tar.xz` so CRAN can build offline; developer clones
  without vendoring build online against the in-repo `rsxcore` path
  dependency.

## Related packages

Downstream R analyses on the same marker outputs are developed in
companion packages (ChromSex, SexPCR) that consume RADSex-compatible
TSVs. Use the CLI, Python (`pyrsx`), or this package to produce those
tables.

## Further reading

- Project site: <https://rsx.rgoswami.me>
- Source and issues: <https://github.com/HaoZeke/rsx-rs>
- Contributing: repository root `CONTRIBUTING.md` and package
  `README.md`
