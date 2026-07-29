# Create a marker table handle

Create a marker table handle

## Usage

``` r
marker_table(x)
```

## Arguments

- x:

  A path to an rsx marker depth table (TSV), or a data frame in the rsx
  marker-table layout.

## Value

An object of class `marker_table`.

## Examples

``` r
tmp <- tempfile(fileext = ".tsv")
writeLines(c("#Number of markers: 0", "id\tsequence"), tmp)
mt <- marker_table(tmp)
print(mt)
#> <marker_table>
#>   path: /tmp/RtmpBAjeAE/file27032fc18ae8.tsv 
```
