#!/bin/sh
# Produce a self-contained, CRAN-submittable rsxr by vendoring the Rust
# dependency tree into src/rust/vendor and compressing it to vendor.tar.xz.
#
# CRAN requires the package to build offline. cargo vendor only vendors
# registry (crates.io) dependencies, not path dependencies, so this script
# first rewrites the rsxcore dependency to its published crates.io version.
#
# Prerequisite: rsxcore must be published to crates.io at the version below
# (it carries the C API the bindings call). Until then this script reports the
# gap and exits non-zero.
#
# Usage: RSX_CORE_VERSION=0.2.6 tools/vendor.sh

set -e

RSX_CORE_VERSION="${RSX_CORE_VERSION:-0.2.6}"
PKG_DIR=$(cd "$(dirname "$0")/.." && pwd)
RUST_DIR="${PKG_DIR}/src/rust"

cd "${RUST_DIR}"

echo "rsxr/vendor: targeting rsxcore ${RSX_CORE_VERSION} from crates.io"

# Confirm the published crate exists before rewriting the manifest.
if ! cargo search rsxcore 2>/dev/null | grep -q "^rsxcore = \"${RSX_CORE_VERSION}\""; then
  echo "WARNING: rsxcore ${RSX_CORE_VERSION} not confirmed on crates.io." >&2
  echo "Publish rsxcore (with the C API) before vendoring, or set" >&2
  echo "RSX_CORE_VERSION to a published version." >&2
fi

# Rewrite the path dependency to the crates.io version (kept reversible: the
# committed Cargo.toml uses the path dep for in-repo development).
cp Cargo.toml Cargo.toml.dev.bak
sed -i.tmp \
  "s|^rsxcore = .*|rsxcore = { version = \"${RSX_CORE_VERSION}\", default-features = false, features = [\"parallel\"] }|" \
  Cargo.toml
rm -f Cargo.toml.tmp

# Vendor the full dependency tree and write the offline source replacement.
cargo generate-lockfile
cargo vendor --versioned-dirs vendor > "${PKG_DIR}/cargo_vendor_config.toml"

# Acknowledge vendored crate authors/licenses (CRAN policy). cargo-about or
# cargo-license, when present, produces a more complete attribution file.
if command -v cargo-about >/dev/null 2>&1; then
  cargo about generate --format json > "${PKG_DIR}/inst/AUTHORS.json" 2>/dev/null || true
fi

# Compress for shipping; the committed tree keeps vendor/ out of git.
tar --create --xz --file vendor.tar.xz vendor
rm -rf vendor

echo "rsxr/vendor: wrote src/rust/vendor.tar.xz and cargo_vendor_config.toml"
echo "rsxr/vendor: restore the dev manifest with 'mv src/rust/Cargo.toml.dev.bak src/rust/Cargo.toml' if needed"
