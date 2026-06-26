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
