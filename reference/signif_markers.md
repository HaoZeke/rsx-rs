# Significant sex-linked markers

Significant sex-linked markers

## Usage

``` r
signif_markers(x, ...)

# S3 method for class 'marker_table'
signif_markers(x, popmap, ...)
```

## Arguments

- x:

  A
  [marker_table](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md).

- ...:

  Passed to
  [`rsx_signif()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_signif.md).

- popmap:

  Population map path.

## Value

A tibble of significant markers.

## Examples

``` r
if (FALSE) { # \dontrun{
mt <- marker_table("markers.tsv")
signif_markers(mt, popmap = "popmap.tsv", test = "fisher")
} # }
```
