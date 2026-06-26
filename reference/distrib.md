# Marker distribution across groups

Marker distribution across groups

## Usage

``` r
distrib(x, ...)

# S3 method for class 'marker_table'
distrib(x, popmap, ...)
```

## Arguments

- x:

  A
  [marker_table](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md).

- ...:

  Passed to
  [`rsx_distrib()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_distrib.md).

- popmap:

  Population map path.

## Value

A tibble of the marker distribution.

## Examples

``` r
# \donttest{
mt <- marker_table("markers.tsv")
#> Error: marker_table: file does not exist: markers.tsv
distrib(mt, popmap = "popmap.tsv", group1 = "M", group2 = "F")
#> Error: object 'mt' not found
# }
```
