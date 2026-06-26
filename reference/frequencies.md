# Per-marker allele frequencies

Per-marker allele frequencies

## Usage

``` r
frequencies(x, ...)

# S3 method for class 'marker_table'
frequencies(x, ...)
```

## Arguments

- x:

  A
  [marker_table](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md).

- ...:

  Passed to
  [`rsx_freq()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_freq.md).

## Value

A tibble of frequencies.

## Examples

``` r
# \donttest{
mt <- marker_table("markers.tsv")
#> Error: marker_table: file does not exist: markers.tsv
frequencies(mt, min_depth = 5L)
#> Error: object 'mt' not found
# }
```
