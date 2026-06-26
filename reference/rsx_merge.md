# Merge marker depth tables (bounded memory)

Merge marker depth tables (bounded memory)

## Usage

``` r
rsx_merge(input_files, output_file, buffer_size = 0L, output_parquet = FALSE)
```

## Arguments

- input_files:

  Character vector of marker table paths.

- output_file:

  Output path.

- buffer_size:

  External-sort buffer size; 0 uses the default.

- output_parquet:

  Write Parquet instead of TSV.

## Value

The output path, invisibly.
