# rsxr

<!-- badges: start -->
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

R bindings for [rsx](https://github.com/HaoZeke/rsx-rs), a high-performance
streaming toolkit for RAD-seq sex determination.

rsxr talks to rsx through R's native C interface: the package's C glue calls
the rsx C API directly and links against the rsxcore static library. There is
no subprocess and no intermediate serialization beyond the marker TSVs the
commands already stream.

## Requirements

- A Rust toolchain (`cargo`, `rustc`); see <https://rustup.rs>.
- The rsxcore sources. Inside the rsx-rs repository these are found
  automatically at `../rsxcore`; otherwise set the `RSX_CORE_DIR` environment
  variable before installing.

## Install

```r
# from inside the rsx-rs repository
install.packages("rsx-r", repos = NULL, type = "source")
```

or with pixi (reproducible R + cargo toolchain):

```bash
cd rsx-r
pixi run install
```

## Usage

### Low-level (mirrors the CLI)

```r
library(rsxr)

rsx_process("reads/", "markers.tsv", threads = 8, min_depth = 5)
rsx_signif("markers.tsv", "popmap.tsv", "signif.tsv",
           test = "fisher", correction = "fdr", bayes = TRUE)
rsx_triage("markers.tsv", "popmap.tsv", "triage.tsv", min_depth = 10)
```

### High-level (tibbles)

```r
mt <- marker_table("markers.tsv")

triaged <- triage(mt, popmap = "popmap.tsv", min_depth = 10)   # tibble
sig     <- signif_markers(mt, popmap = "popmap.tsv", test = "fisher")
dist    <- distrib(mt, popmap = "popmap.tsv", group1 = "M", group2 = "F")
freqs   <- frequencies(mt, min_depth = 5)
```

## Ecosystem

rsxr is the R entry point to the rsx engine; the analysis-layer R packages
[ChromSex](https://github.com/RuhiRG/ChromSex) and
[SexPCR](https://github.com/RuhiRG/SexPCR) build on the same marker outputs.

## License

GPL-3.0-or-later. See `LICENSE`.
