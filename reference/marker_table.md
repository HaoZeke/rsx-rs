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
if (FALSE) { # \dontrun{
mt <- marker_table("markers.tsv")
triage(mt, popmap = "popmap.tsv", min_depth = 10)
} # }
```
