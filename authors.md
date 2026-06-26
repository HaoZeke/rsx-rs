# Authors and Citation

## Authors

- **Rohit Goswami**. Author, maintainer.
  [](https://orcid.org/0000-0002-2393-8056)

- **Ruhila Goswami**. Author. [](https://orcid.org/0000-0002-5443-9356)

## Citation

Source:
[`inst/CITATION`](https://github.com/HaoZeke/rsx-rs/blob/main/inst/CITATION)

Goswami R, Goswami R (2026). rsxr: R Bindings for the rsx RAD-seq
Sex-Determination Toolkit. https://github.com/HaoZeke/rsx-rs

    @Manual{,
      title = {rsxr: R Bindings for the rsx RAD-seq Sex-Determination Toolkit},
      author = {Rohit Goswami and Ruhila Goswami},
      year = {2026},
      note = {R package version 0.1.0},
      url = {https://github.com/HaoZeke/rsx-rs},
    }

## Additional details

    rsxr links against the rsxcore Rust library and, in a CRAN build, bundles its
    Rust dependency tree under src/rust/vendor (compressed as vendor.tar.xz).

    The bundled crates are third-party Rust packages distributed under permissive
    licenses (predominantly MIT or Apache-2.0). A machine-readable manifest of the
    exact vendored crates, versions, and licenses is produced by tools/vendor.sh
    (via cargo-about, written to inst/AUTHORS.json) at packaging time and ships in
    the source tarball.

    rsxcore and these bindings are released under GPL-3.0-or-later.

    Upstream:
    - rsxcore: https://github.com/HaoZeke/rsx-rs (Rohit Goswami, Ruhila Goswami)
    - needletail, flate2, rayon, and the rest of the dependency tree: see the
      vendored Cargo.lock and inst/AUTHORS.json.
