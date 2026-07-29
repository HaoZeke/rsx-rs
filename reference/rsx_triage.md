# Bayesian sex-linkage triage of markers

Bayesian sex-linkage triage of markers

## Usage

``` r
rsx_triage(
  table_path,
  popmap_path,
  output_file,
  min_depth = 10L,
  signif_threshold = 0.05,
  posterior_threshold = 0.9,
  bayes_factor_threshold = 10,
  prior_probability = 0.01,
  linked_probability = 0.9,
  null_prevalence = 0.5,
  group1_linked_weight = 0.5,
  group1 = "M",
  group2 = "F"
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

- posterior_threshold:

  Posterior probability cutoff.

- bayes_factor_threshold:

  Bayes factor cutoff.

- prior_probability:

  Prior probability that a marker is sex-linked.

- linked_probability:

  Expected marker prevalence in the linked group.

- null_prevalence:

  Expected marker prevalence under the null model.

- group1_linked_weight:

  Mixture weight for the group-1-linked direction.

- group1, group2:

  Group labels to contrast (e.g. "M", "F").

## Value

The output path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
rsx_triage("markers.tsv", "popmap.tsv", tempfile(fileext = ".tsv"),
           min_depth = 10L)
} # }
```
