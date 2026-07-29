# rsxr: R bindings for the rsx RAD-seq sex-determination toolkit

rsxr binds the rsx C API directly through R's native C interface. The
low-level functions
([`rsx_process()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_process.md),
[`rsx_freq()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_freq.md),
[`rsx_distrib()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_distrib.md),
[`rsx_signif()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_signif.md),
[`rsx_triage()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_triage.md),
[`rsx_depth()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_depth.md),
[`rsx_merge()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_merge.md),
[`rsx_pca()`](https://haozeke.github.io/rsx-rs/rsxr/reference/rsx_pca.md))
mirror the rsx CLI one to one. The high-level
[`marker_table()`](https://haozeke.github.io/rsx-rs/rsxr/reference/marker_table.md)
workflow returns tibbles.

## See also

Useful links:

- <https://github.com/HaoZeke/rsx-rs>

- <https://rsx.rgoswami.me>

- Report bugs at <https://github.com/HaoZeke/rsx-rs/issues>

## Author

**Maintainer**: Rohit Goswami <rgoswami@ieee.org>
([ORCID](https://orcid.org/0000-0002-2393-8056))

Authors:

- Rohit Goswami <rgoswami@ieee.org>
  ([ORCID](https://orcid.org/0000-0002-2393-8056))

- Ruhila Goswami ([ORCID](https://orcid.org/0000-0002-5443-9356))
