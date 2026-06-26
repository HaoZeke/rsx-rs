# Read an rsx core TSV output

rsx commands prefix their TSV output with a `#Number of markers` comment
line; this reads past it and returns a tibble.

## Usage

``` r
rsxr_read_core_tsv(path)
```

## Arguments

- path:

  TSV path.

## Value

A tibble.
