<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/branding/rsx-logo-dark.svg">
    <img src="assets/branding/rsx-logo.svg" alt="rsx" width="460">
  </picture>
</p>

`rsx` builds a marker depth table from demultiplexed RAD-seq reads and tests markers for association with two groups named in a population map. The commands follow [RADSex](https://github.com/RomainFeron/RADSex) (`process`, `distrib`, `signif`, `freq`, `depth`, `map`, `subset`) and add `merge`, `pca`, and `triage`. The library crate is `rsxcore`. The binary is `rsx` (`rsx-cli`). Python bindings are the `pyrsx` package.

[![CI](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/HaoZeke/rsx-rs/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/HaoZeke/rsx-rs/graph/badge.svg)](https://app.codecov.io/gh/HaoZeke/rsx-rs)
[![Documentation](https://img.shields.io/badge/docs-rsx.rgoswami.me-blue)](https://rsx.rgoswami.me)
[![Crates.io](https://img.shields.io/crates/v/rsx-cli?label=crates.io)](https://crates.io/crates/rsx-cli)
[![PyPI](https://img.shields.io/pypi/v/pyrsx)](https://pypi.org/project/pyrsx/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Runtime comparisons with C++ RADSex are in the paper cited below. Command output matches RADSex when the same groups are given.

## Install

Rust 1.85 or newer (edition 2024). The published binary is 0.2.10.

### GitHub release

Pre-built binaries for Linux (x86_64 and aarch64), macOS (x86_64 and arm64), and Windows (the Windows build omits `map`):

```bash
# https://github.com/HaoZeke/rsx-rs/releases
curl --proto '=https' --tlsv1.2 -LsSf https://github.com/HaoZeke/rsx-rs/releases/download/v0.2.10/rsx-cli-installer.sh | sh
```

### From source

```bash
git clone https://github.com/HaoZeke/rsx-rs.git
cd rsx-rs
cargo build --release
# binary: target/release/rsx
```

`map` is enabled by default on the CLI. CUDA chi-square batches and parquet output are off unless requested:

```bash
cargo build -p rsx-cli --release --features cuda
cargo build -p rsx-cli --release --features parquet-io
```

A build that does not include the feature returns an error if you pass `--backend cuda` or `--output-parquet`. It does not panic.

### pixi

```bash
pixi run build
# or: pixi run -e dev build-portable
```

### Python

```bash
pip install pyrsx
```

The [Python README](rsx-python/README.md) documents `MarkerTable`, `TriageResult`, narwhals, and plotting. R bindings are in [rsx-r/README.md](rsx-r/README.md).

## Worked example

The commands below use the committed panel in `benchmarks/data/small`: ten gzipped FASTQ files (`ind1` through `ind10`), `popmap.tsv` (`ind1`–`ind5` are `M`, `ind6`–`ind10` are `F`), and `genome.fa`. Individual names are the filenames with the sequence extension removed, sorted lexicographically, so the depth columns are `ind1`, `ind10`, `ind2`, …, `ind9`.

```bash
rsx process -i benchmarks/data/small/reads -o markers.tsv -T 2 -d 5
```

`process` keeps a marker when at least one individual has depth `>= 5`. This panel writes 469 markers:

```text
#Number of markers : 469
id	sequence	ind1	ind10	ind2	ind3	ind4	ind5	ind6	ind7	ind8	ind9
0	AAAAATTCGAGATACCTTTTACCAGGCAGAGACAAAAGTTCGGATGACAGATATTCCAGTGACAAAAGAGACATCACCAAGTTAAAAATGGACTGTTATC	11	0	0	6	4	0	0	0	0	5
1	AAACATGAAGTTAACGTCTGCTCATGACACGCACGACTCGACAACTGGGGGTCGATGCACCAACGCCTCGGAGTGTTAAACGGTCCTACCAAAAAGGTGC	0	8	0	0	0	0	5	10	3	8
```

Frequency histogram (`freq` does not read the popmap). Frequency 0 is omitted:

```bash
rsx freq -t markers.tsv -o freq.tsv -d 5
```

```text
#source:rsx-freq;min_depth:5
Frequency	Count
1	51
2	97
3	131
4	115
5	52
6	20
7	2
8	1
9	0
10	0
```

Distribution of presence counts for groups `M` and `F`. The comment records the Bonferroni threshold `0.05 / 469`:

```bash
rsx distrib -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -o distrib.tsv -G M,F -d 5
```

```text
#source:rsx-distrib;min_depth:5;signif_threshold:0.00010661;correction:bonferroni;n_markers:469;n_tests:469
M	F	Markers	P	CorrectedP	Signif	Bias
0	1	13	1	1	False	-0.2
0	2	18	0.429195	1	False	-0.4
0	3	18	0.167546	1	False	-0.6
0	4	26	0.0528075	1	False	-0.8
0	5	10	0.011412	1	False	-1
1	0	38	1	1	False	0.2
```

The full table has 35 cells (every pair with a defined test, including cells whose marker count is 0). None of them pass Bonferroni on this panel.

`signif` with the default correction (`bonferroni`) and `--bayes` writes the header and no data rows:

```bash
rsx signif -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -o signif.tsv -G M,F -d 5 --bayes
```

```text
#source:rsx-signif;min_depth:5;signif_threshold:0.05;correction:bonferroni;test:chisq;n_markers:469
id	sequence	ind1	ind10	ind2	ind3	ind4	ind5	ind6	ind7	ind8	ind9	Bayes_Factor	Posterior_SexLinked
```

`--correction none` keeps markers whose raw p-value is below 0.05. That run writes 10 data rows. The first is:

```bash
rsx signif -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -o signif-none.tsv -G M,F -d 5 --bayes --correction none
```

```text
#source:rsx-signif;min_depth:5;signif_threshold:0.05;correction:none;test:chisq;n_markers:469
id	sequence	ind1	ind10	ind2	ind3	ind4	ind5	ind6	ind7	ind8	ind9	Bayes_Factor	Posterior_SexLinked
144	CAGGGATTGAACCACGGGCCTCCTCCGGAAAACGCATGGCAACGGTTCGCCATCTAGGTTCACGCCTGAGCTGTCTGTGACGCAGATTGCGGCTCGGTCT	0	8	0	0	0	0	6	8	6	6	77.0000	0.6433
```

Per-individual depth summary. `Retained` is the number of markers present in at least 75% of individuals (`--min-frequency`, default 0.75):

```bash
rsx depth -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -o depth.tsv
```

```text
Sample	Group	Reads	Markers	Retained	Min_depth	Max_depth	Median_depth	Average_depth
ind1	M	1451	299	207	0	17	4	4
ind10	F	1507	290	207	0	18	3	3
ind2	M	1468	296	207	0	25	3	3
ind3	M	1453	290	207	0	14	4	3
ind4	M	1438	303	207	0	13	4	3
ind5	M	1507	296	207	0	18	4	4
ind6	F	1489	293	207	0	22	3	3
ind7	F	1433	291	207	0	14	2	3
ind8	F	1476	282	207	0	16	3	3
ind9	F	1485	296	207	0	18	4	3
```

`triage` writes one row per retained marker (41 lines: comment, header, 39 markers). The first data row on this panel is class `bayes_factor_only`:

```bash
rsx triage -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -o triage.tsv -G M,F -d 5
```

```text
#source:rsx-triage;min_depth:5;signif_threshold:0.05;posterior_threshold:0.9;bayes_factor_threshold:10;prior_probability:0.01;linked_probability:0.9;n_markers:469
id	sequence	Group1	Group1_Present	Group1_Total	Group1_Penetrance	Group2	Group2_Present	Group2_Total	Group2_Penetrance	Bias_Direction	Bias	P	CorrectedP	Bayes_Factor	Posterior_SexLinked	Strict_Call	Posterior_Call	Bayes_Factor_Call	Candidate_Class
1	AAACATGAAGTTAACGTCTGCTCATGACACGCACGACTCGACAACTGGGGGTCGATGCACCAACGCCTCGGAGTGTTAAACGGTCCTACCAAAAAGGTGC	M	0	5	0	F	4	5	0.8	F-biased	-0.8	0.0528075	1	12.8333	0.166919	false	false	true	bayes_factor_only
```

Sample-space PCA of the depth matrix, two components:

```bash
rsx pca -t markers.tsv -o pca/ -d 5 -r 2
```

`pca/summary.txt`:

```text
Streaming PCA of depth matrix
Markers: 469
Individuals: 10
Components: 2
Total variance: 55294.83

Top components:
  PC1: 0.1890 variance (18.9% cumulative)
  PC2: 0.1187 variance (30.8% cumulative)
```

`pca/eigenvalues.tsv`:

```text
component	eigenvalue	variance_fraction	cumulative
PC1	10451.027234	0.189005	0.189005
PC2	6560.780864	0.118651	0.307656
```

`pca/loadings.tsv` and `pca/sample_scores.tsv` have the same contents (`individual`, `PC1`, `PC2`, …).

Align markers to the bundled genome (`map` needs the `map` feature, which the default CLI build has):

```bash
rsx map -t markers.tsv -p benchmarks/data/small/popmap.tsv \
  -g benchmarks/data/small/genome.fa -o map.tsv -G M,F -d 5
```

```text
#source:rsx-map;min_depth:5;min_qual:20;min_freq:0.1;signif_threshold:0.00010661;bonferroni:true;n_markers:469
Contig	Position	Length	Marker_id	Bias	P	CorrectedP	Signif
chr3	2307	10000	0	0.2	1	1	False
chr1	8605	10000	1	-0.8	0.0528075	1	False
```

An unknown test name exits 1 and prints the accepted values. It does not panic:

```text
$ rsx signif -t markers.tsv -p benchmarks/data/small/popmap.tsv -o out.tsv -G M,F --test anova
error: signif: test 'anova' is not recognized. Accepted values: chisq (aliases: chi-squared, chi2), fisher (alias: fisher-exact), gtest (aliases: g-test, g)
```

`rsx --help` lists commands. `rsx <command> --help` lists that command's arguments and the columns it writes. `-h` is the short summary. `-V`, `--version`, and `-v` print the version and exit 0.

## Inputs

### Reads (`process -i`)

A directory of demultiplexed FASTA or FASTQ, optionally gzipped. Accepted extensions: `.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`, `.fa`, `.fa.gz`, `.fasta`, `.fasta.gz`, `.fna`, `.fna.gz`. One file is one individual. `--threads` / `-T` must be at least 1. `--kmer-dedup` / `-k`, when set, must be at least 1; omit it to disable k-mer deduplication.

### Population map

One `individual<TAB>group` row per line. Blank lines are ignored. A line that does not have two columns is skipped and a warning is written to stderr. A popmap that names zero individuals is accepted: `depth` then leaves `Group` blank. `--groups` / `-G` takes exactly two non-empty names, for example `--groups M,F`. If the popmap contains exactly two groups, `distrib`, `signif`, `triage`, `subset`, and `map` can omit `--groups`.

### Marker depth table

Written by `process` and read by every later command:

```text
#Number of markers : N
id	sequence	<individual>	<individual>	...
<id>	<sequence>	<depth>	<depth>	...
```

The `#Number of markers` line is optional on input. The header must start with the columns `id` and `sequence`, then one column per individual. Depths are integers. `merge` saturates a summed depth at 65535.

## Outputs

| Command | File | Contents |
| --- | --- | --- |
| `process` | TSV | Marker depth table. Comment `#Number of markers : N`, then `id`, `sequence`, one column per individual. |
| `freq` | TSV | `#source:rsx-freq;min_depth:N`, then `Frequency` and `Count`. Frequency 0 is omitted. |
| `distrib` | TSV | Comment, then `<group1>`, `<group2>`, `Markers`, `P`, `CorrectedP`, `Signif`, `Bias`. `--bayes` appends `Bayes_Factor` and `Posterior_SexLinked` (four decimals). |
| `signif` | TSV or FASTA | Comment, then the input header and the rows that pass the corrected threshold. `--bayes` appends the same two columns. FASTA headers are `>id_<group1>:<count>_<group2>:<count>_p:<p>_pcorr:<pcorr>_mindepth:<d>`. `--output-fasta` with `--correction fdr` is an error. |
| `triage` | TSV | One row per retained marker: group counts, penetrance, bias, p, corrected p, Bayes factor, posterior, three boolean calls, and `Candidate_Class` (`strict+posterior`, `strict_only`, `posterior_only`, `bayes_factor_only`, `exploratory`). The strict call uses `signif_threshold / n_markers`. |
| `depth` | TSV | `Sample`, `Group`, `Reads`, `Markers`, `Retained`, `Min_depth`, `Max_depth`, `Median_depth`, `Average_depth`. |
| `subset` | TSV or FASTA | Depth table with a `#source:rsx-subset` comment, or FASTA in the `signif` header shape. Default maxima of 9999 leave those bounds unset for typical cohorts. |
| `merge` | TSV or parquet | Depth tables concatenated in the order given. Equal sequences have depths added with saturation at 65535. Progress is written to stderr and is not part of the table. `--output-parquet` needs the `parquet-io` feature. |
| `pca` | directory | `eigenvalues.tsv`, `loadings.tsv`, `sample_scores.tsv` (same rows as `loadings.tsv`), and `summary.txt`. |
| `map` | TSV | `Contig`, `Position`, `Length`, `Marker_id`, `Bias`, `P`, `CorrectedP`, `Signif`. One row per alignment with mapping quality at least `--min-quality` (default 20). Minimum allele frequency default is 0.1. |

Tests: `chisq` (aliases `chi-squared`, `chi2`; default), `fisher` (`fisher-exact`), `gtest` (`g-test`, `g`). Corrections: `bonferroni` (`bonf`; default), `fdr` (`bh`, `benjamini-hochberg`), `none` (`disabled`). Backend: `cpu` (default) or `cuda` (`gpu`).

`depth` chooses memory or an external sort from the file size (`--streaming-mode auto|memory|streaming`, `--streaming-threshold-bytes`). The median is the same in both modes. FDR `signif` stores one p-value per marker for Benjamini-Hochberg, then streams the table again. The other default paths stream the depth table.

## Run profiles

A profile is TOML with `schema_version = 1` and one `[run]` command. Flags on the command line override the profile. The hydrated profile is the values the calculation receives, including Bayesian prior shapes and the depth streaming policy.

```bash
rsx --profile examples/profiles/triage.toml \
  --write-hydrated-profile triage.hydrated.toml \
  --reproducibility-archive triage-reproducibility.zip
```

The reproducibility zip is created before input validation. A failed run updates that zip when profile resolution succeeded. The zip holds the input and hydrated TOML, the executable, the lockfile, build and run manifests, JSON Schema, a CycloneDX software bill of materials, citation and license files, and checksums. It does not hold analysis results or input datasets.

## Errors

Failures go to stderr as `error: ...` and exit 1. The message names the command, the path when a file could not be opened, and the accepted values when a flag is wrong. On this panel:

```text
error: process: n_threads must be at least 1 (CLI: --threads / -T)
error: process: kmer_dedup must be at least 1 when set (CLI: --kmer-dedup). Omit it to disable deduplication.
error: cannot open markers table /tmp/no-such-markers.tsv: No such file or directory (os error 2). Expected a TSV with an optional '#Number of markers : N' line and a header row 'id\tsequence\t<individual>...'
error: --groups requires exactly two non-empty group names separated by a comma, for example --groups M,F; got 1 (M)
error: backend cuda was requested, but this build does not include the cuda feature. Rebuild with `cargo build -p rsx-cli --release --features cuda`, or pass --backend cpu
```

When `--reproducibility-archive` is set, these errors are recorded on the archive after the message is printed.

## Documentation

Command reference, tutorials, and the memory table: https://rsx.rgoswami.me

Library docs for the current tree: `cargo doc --open -p rsxcore --no-deps`.

Paper: Goswami R, Goswami R. rsx: a high-performance streaming toolkit for RAD-seq sex determination. *BMC Bioinformatics* (2026). https://doi.org/10.1186/s12859-026-06628-4

Preprint: https://arxiv.org/abs/2606.06434

Software archive: https://doi.org/10.5281/zenodo.20531538

Downstream of the marker table, [ChromSex](https://github.com/RuhiRG/ChromSex) computes sex-balanced Hudson / Weir and Cockerham FST, sex-specific heterozygosity, and genomic content around sex loci. [SexPCR](https://github.com/RuhiRG/SexPCR) screens sdY PCR primers against dissection sex.

## Citation

Cite the BMC Bioinformatics paper.

> Goswami R, Goswami R. rsx: a high-performance streaming toolkit for RAD-seq
> sex determination. BMC Bioinformatics (2026).
> https://doi.org/10.1186/s12859-026-06628-4

`CITATION.cff` is the machine-readable record (journal DOI, arXiv preprint, RADSex, and the Zenodo software DOI). GitHub's "Cite this repository" button reads it.

RADSex: Feron et al., Molecular Ecology Resources 2021. https://doi.org/10.1111/1755-0998.13360

Benchmark data and figures: https://doi.org/10.5281/zenodo.20531539

### BibTeX

```bibtex
@article{goswamiRsxHighperformanceStreaming2026,
  title         = {rsx: a high-performance streaming toolkit for {RAD-seq} sex determination},
  author        = {Goswami, Rohit and Goswami, Ruhila},
  year          = {2026},
  month         = sep,
  journal       = {BMC Bioinformatics},
  publisher     = {BioMed Central},
  doi           = {10.1186/s12859-026-06628-4},
  url           = {https://doi.org/10.1186/s12859-026-06628-4},
  eprint        = {2606.06434},
  archivePrefix = {arXiv},
  primaryClass  = {q-bio.GN}
}
```

For the software archive (v0.2.10) use `CITATION.cff` or the Zenodo record. The reproducibility archive has its own Zenodo BibTeX.

## Contributing

See `CONTRIBUTING.md`. Development uses pixi, conventional commits, and:

```bash
cargo fmt --all
cargo clippy --workspace -- -D warnings
cargo test --workspace
```

## License

GPL-3.0-or-later. See `LICENSE`.
