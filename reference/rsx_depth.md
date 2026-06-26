# Marker depth per individual

Marker depth per individual

## Usage

``` r
rsx_depth(
  table_path,
  popmap_path,
  output_file,
  min_frequency = 0.5,
  streaming = FALSE
)
```

## Arguments

- table_path:

  Marker depth table (TSV).

- popmap_path:

  Population map (TSV).

- output_file:

  Output path (TSV).

- min_frequency:

  Minimum group presence frequency.

- streaming:

  Use the streaming (bounded-memory) path.

## Value

The output path, invisibly.

## Examples

``` r
# \donttest{
rsx_depth("markers.tsv", "popmap.tsv", tempfile(fileext = ".tsv"))
#> Error in rsx_depth("markers.tsv", "popmap.tsv", tempfile(fileext = ".tsv")): rsx: depth failed: No such file or directory (os error 2) (status 2)
# }
```
