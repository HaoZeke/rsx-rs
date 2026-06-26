# Streaming PCA of the marker matrix

Streaming PCA of the marker matrix

## Usage

``` r
rsx_pca(table_path, output_dir, min_depth = 1L, n_components = 0L)
```

## Arguments

- table_path:

  Marker depth table (TSV).

- output_dir:

  Directory for PCA outputs.

- min_depth:

  Minimum depth threshold.

- n_components:

  Number of components; 0 uses the default.

## Value

The output directory, invisibly.
