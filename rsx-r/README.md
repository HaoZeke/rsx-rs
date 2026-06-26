# rsxr

<!-- badges: start -->
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/HaoZeke/rsx-rs/actions/workflows/rsxr-R-CMD-check.yaml/badge.svg)](https://github.com/HaoZeke/rsx-rs/actions/workflows/rsxr-R-CMD-check.yaml)
[![pkgcheck](https://github.com/HaoZeke/rsx-rs/actions/workflows/rsxr-pkgcheck.yaml/badge.svg)](https://github.com/HaoZeke/rsx-rs/actions/workflows/rsxr-pkgcheck.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

R bindings for [rsx](https://github.com/HaoZeke/rsx-rs), a high-performance
streaming toolkit for RAD-seq sex determination.

## Statement of need

Sex-determination studies from RAD-seq data need marker depth tables, group
contrasts, and frequentist or Bayesian triage. The **rsx** engine implements
those steps with bounded-memory streaming and a C API. **rsxr** exposes that
API in-process from R (no subprocess), so analysts can run `process`, `freq`,
`distrib`, `signif`, `triage`, `depth`, `merge`, and `pca` and consume results
as tibbles without leaving R. It targets R-first workflows that already use
popmaps and RADSex-compatible marker TSVs.

rsxr talks to rsx through R's native C interface: the package's C glue calls
the rsx C API directly and links against the rsxcore static library. There is
no intermediate serialization beyond the marker TSVs the commands already stream.

## Requirements

- A Rust toolchain (`cargo`, `rustc`); see <https://rustup.rs>. Minimum supported
  Rust version is recorded in `DESCRIPTION` as `Config/rsxr/MSRV`.
- The **rsxcore** sources as a sibling of this package (`../rsxcore` inside the
  rsx-rs monorepo). Otherwise set `RSX_CORE_DIR` before installing.
- On Linux, development headers for zlib, bzip2, and xz if the linker cannot
  find them.

## Install

From GitHub (package lives in the `rsx-r` subdirectory):

```r
# pak
pak::pak("HaoZeke/rsx-rs/rsx-r")

# remotes
remotes::install_github("HaoZeke/rsx-rs", subdir = "rsx-r")
```

From a local clone of [rsx-rs](https://github.com/HaoZeke/rsx-rs):

```bash
cd rsx-r
R CMD INSTALL .
# or with the package pixi environment (R + cargo pinned):
pixi run install
```

## Usage

### Low-level (mirrors the CLI)

```r
library(rsxr)

rsx_version()

rsx_process("reads/", "markers.tsv", threads = 8L, min_depth = 5L)
rsx_signif("markers.tsv", "popmap.tsv", "signif.tsv",
           test = "fisher", correction = "fdr", bayes = TRUE)
rsx_triage("markers.tsv", "popmap.tsv", "triage.tsv", min_depth = 10L)
```

### High-level (tibbles)

```r
mt <- marker_table("markers.tsv")

triaged <- triage(mt, popmap = "popmap.tsv", min_depth = 10L)
sig     <- signif_markers(mt, popmap = "popmap.tsv", test = "fisher")
dist    <- distrib(mt, popmap = "popmap.tsv", group1 = "M", group2 = "F")
freqs   <- frequencies(mt, min_depth = 5L)
```

See `vignette("rsxr", package = "rsxr")` after install.

## Ecosystem

rsxr is the R entry point to the rsx engine. Companion R packages for
sex-balanced FST / heterozygosity (ChromSex) and sdY PCR assay audit (SexPCR)
consume the same marker outputs; those repositories are not yet public.

## Contributing and community guidelines

Contributions are welcome. Please open issues and pull requests on the
[rsx-rs repository](https://github.com/HaoZeke/rsx-rs).

- Contribution guidelines: repository root
  [`CONTRIBUTING.md`](https://github.com/HaoZeke/rsx-rs/blob/main/CONTRIBUTING.md)
- Code of conduct: repository root
  [`CODE_OF_CONDUCT.md`](https://github.com/HaoZeke/rsx-rs/blob/main/CODE_OF_CONDUCT.md)
  (if present) or the GitHub Community Guidelines

For the R package specifically:

```bash
cd rsx-r
pixi install
pixi run document   # roxygen2
pixi run test
pixi run check      # rcmdcheck --no-manual
```

## License

GPL-3. See `LICENSE` and `LICENSE.md`.
