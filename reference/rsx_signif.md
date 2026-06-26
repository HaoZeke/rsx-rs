# Extract significant sex-linked markers

Extract significant sex-linked markers

## Usage

``` r
rsx_signif(
  table_path,
  popmap_path,
  output_file,
  min_depth = 1L,
  signif_threshold = 0.05,
  group1 = "",
  group2 = "",
  correction = "bonferroni",
  test = "chisq",
  output_fasta = FALSE,
  bayes = FALSE
)
```

## Arguments

- table_path:

  Marker depth table (TSV).

- popmap_path:

  Population map (TSV).

- output_file:

  Output path (TSV).

- min_depth:

  Minimum depth threshold.

- signif_threshold:

  Significance threshold.

- group1, group2:

  Group labels to contrast (e.g. "M", "F").

- correction:

  Multiple-testing correction: "bonferroni", "fdr", "none".

- test:

  Statistical test: "chisq", "fisher", "gtest".

- output_fasta:

  Also write a FASTA of significant markers.

- bayes:

  Also emit Bayesian columns.

## Value

The output path, invisibly.
