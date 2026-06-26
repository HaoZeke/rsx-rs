# Marker distribution across groups

Marker distribution across groups

## Usage

``` r
rsx_distrib(
  table_path,
  popmap_path,
  output_file,
  min_depth = 1L,
  signif_threshold = 0.05,
  group1 = "",
  group2 = "",
  correction = "bonferroni",
  test = "chisq",
  output_bayes = FALSE
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

- output_bayes:

  Also emit Bayesian columns.

## Value

The output path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
rsx_distrib("markers.tsv", "popmap.tsv", tempfile(fileext = ".tsv"),
            group1 = "M", group2 = "F")
} # }
```
