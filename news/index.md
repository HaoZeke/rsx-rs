# Changelog

## rsxr 0.1.0

#### Added

- Raw C bindings to the rsx C API (no subprocess): `rsx_process`,
  `rsx_freq`, `rsx_distrib`, `rsx_signif`, `rsx_triage`, `rsx_depth`,
  `rsx_merge`, `rsx_pca`, and `rsx_version`.
- Ergonomic `marker_table` workflow returning tibbles, with the verbs
  `triage`, `signif_markers`, `distrib`, `depth`, and `frequencies`.
- `configure` + `Makevars` that build the rsxcore static library with
  cargo and link it into the package.
