# Bayesian sex-linkage triage

Bayesian sex-linkage triage

## Usage

``` r
triage(x, ...)

# S3 method for class 'marker_table'
triage(x, popmap, ...)
```

## Arguments

- x:

  A
  [marker_table](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md).

- ...:

  Passed to
  [`rsx_triage()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_triage.md).

- popmap:

  Population map path.

## Value

A tibble of triaged markers.

## Examples

``` r
# \donttest{
mt <- marker_table("markers.tsv")
#> Error: marker_table: file does not exist: markers.tsv
triage(mt, popmap = "popmap.tsv", min_depth = 10L)
#> Error: object 'mt' not found
# }
```
