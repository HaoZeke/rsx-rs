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
  bayes = FALSE,
  prior_probability = 0.01,
  linked_probability = 0.9,
  null_prevalence = 0.5,
  group1_linked_weight = 0.5,
  bf_group1_alpha = 1,
  bf_group1_beta = 1,
  bf_group2_alpha = 1,
  bf_group2_beta = 1,
  bf_null_alpha = 1,
  bf_null_beta = 1,
  posterior_linked_family = "fixed",
  posterior_linked_probability = linked_probability,
  posterior_linked_alpha = 1,
  posterior_linked_beta = 1,
  posterior_null_family = "fixed",
  posterior_null_probability = null_prevalence,
  posterior_null_alpha = 1,
  posterior_null_beta = 1
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

- prior_probability:

  Prior probability that a marker is sex-linked.

- linked_probability:

  Expected marker prevalence in the linked group.

- null_prevalence:

  Expected marker prevalence under the null model.

- group1_linked_weight:

  Mixture weight for the group-1-linked direction.

- bf_group1_alpha, bf_group1_beta:

  Beta-prior shape parameters for marker prevalence under the
  group-1-linked alternative.

- bf_group2_alpha, bf_group2_beta:

  Beta-prior shape parameters for marker prevalence under the
  group-2-linked alternative.

- bf_null_alpha, bf_null_beta:

  Beta-prior shape parameters for marker prevalence under the null
  model.

- posterior_linked_family, posterior_null_family:

  Posterior prevalence family: `"fixed"` uses its explicit probability;
  `"beta"` integrates over the corresponding Beta shape parameters.

- posterior_linked_probability, posterior_null_probability:

  Fixed posterior prevalence probabilities. Their defaults preserve the
  legacy `linked_probability` and `null_prevalence` behavior.

- posterior_linked_alpha, posterior_linked_beta:

  Beta-prior shape parameters for linked-marker prevalence.

- posterior_null_alpha, posterior_null_beta:

  Beta-prior shape parameters for null-marker prevalence.

## Value

The output path, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
rsx_signif("markers.tsv", "popmap.tsv", tempfile(fileext = ".tsv"),
           test = "fisher", correction = "fdr")
} # }
```
