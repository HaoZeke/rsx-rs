# Compute per-marker allele frequencies

Compute per-marker allele frequencies

## Usage

``` r
rsx_freq(table_path, output_file, min_depth = 1L)
```

## Arguments

- table_path:

  Marker depth table (TSV).

- output_file:

  Output path (TSV).

- min_depth:

  Minimum depth threshold.

## Value

The output path, invisibly.

## Examples

``` r
# \donttest{
rsx_freq("markers.tsv", tempfile(fileext = ".tsv"), min_depth = 5L)
#> Error in rsx_freq("markers.tsv", tempfile(fileext = ".tsv"), min_depth = 5L): rsx: freq failed: No such file or directory (os error 2) (status 2)
# }
```
