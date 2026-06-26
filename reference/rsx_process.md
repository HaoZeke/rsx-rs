# Process demultiplexed reads into a marker depth table

Process demultiplexed reads into a marker depth table

## Usage

``` r
rsx_process(input_dir, output_file, threads = 1L, min_depth = 1L)
```

## Arguments

- input_dir:

  Directory of demultiplexed FASTQ files.

- output_file:

  Path to write the marker depth table (TSV).

- threads:

  Number of worker threads.

- min_depth:

  Minimum depth for a marker to be retained.

## Value

The output path, invisibly.

## Examples

``` r
# \donttest{
rsx_process("reads/", tempfile(fileext = ".tsv"), threads = 2L, min_depth = 5L)
#> Error in rsx_process("reads/", tempfile(fileext = ".tsv"), threads = 2L,     min_depth = 5L): rsx: process failed: No such file or directory (os error 2) (status 2)
# }
```
