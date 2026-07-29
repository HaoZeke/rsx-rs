#!/bin/sh
# Stage monorepo rsxcore into src/rust/rsxcore for package-local cargo builds.
# R CMD build copies only this package, so the path dep cannot point at ../rsxcore.
# Rewrites workspace-inherited fields so the staged crate is a standalone manifest.
set -e
PKG_ROOT="${1:-.}"
DEST="${PKG_ROOT}/src/rust/rsxcore"
CORE_DIR="${RSX_CORE_DIR:-}"

if [ -f "${DEST}/Cargo.toml" ] && [ "${RSXR_FORCE_STAGE:-}" != "1" ]; then
  # Already staged (e.g. CI pre-step). Still ensure standalone manifest.
  :
elif [ -n "${CORE_DIR}" ] && [ -f "${CORE_DIR}/Cargo.toml" ]; then
  :
elif [ -f "${PKG_ROOT}/../rsxcore/Cargo.toml" ]; then
  CORE_DIR="${PKG_ROOT}/../rsxcore"
else
  if [ -f "${DEST}/Cargo.toml" ]; then
    CORE_DIR=""
  else
    echo "rsxcore sources not found (need ${DEST}, RSX_CORE_DIR, or ../rsxcore)" >&2
    exit 1
  fi
fi

if [ -n "${CORE_DIR}" ]; then
  echo "staging rsxcore from ${CORE_DIR} into ${DEST}"
  mkdir -p "${DEST}"
  if command -v rsync >/dev/null 2>&1; then
    rsync -a --delete \
      --exclude target --exclude .git --exclude .pixi \
      "${CORE_DIR}/" "${DEST}/"
  else
    (cd "${CORE_DIR}" && tar cf - \
      --exclude=target --exclude=.git --exclude=.pixi .) |
      (cd "${DEST}" && tar xf -)
  fi
fi

test -f "${DEST}/include/rsx.h"
cp "${DEST}/include/rsx.h" "${PKG_ROOT}/src/rsx.h"

# Read workspace version from monorepo root when available; else default.
WS_VERSION="0.2.6"
if [ -f "${PKG_ROOT}/../Cargo.toml" ]; then
  v=$(sed -n 's/^version = "\([^"]*\)"/\1/p' "${PKG_ROOT}/../Cargo.toml" | head -1)
  [ -n "$v" ] && WS_VERSION="$v"
fi
# Prefer staged package version line if already concrete
if grep -q '^version = "' "${DEST}/Cargo.toml" 2>/dev/null; then
  :
else
  # Rewrite workspace inheritance for a standalone crate manifest.
  tmp="${DEST}/Cargo.toml.stage"
  sed \
    -e "s/^version\\.workspace = true/version = \"${WS_VERSION}\"/" \
    -e "s/^edition\\.workspace = true/edition = \"2024\"/" \
    -e "s/^license\\.workspace = true/license = \"GPL-3.0-or-later\"/" \
    -e "s/^rust-version\\.workspace = true/rust-version = \"1.85\"/" \
    -e "s/^repository\\.workspace = true/repository = \"https:\\/\\/github.com\\/HaoZeke\\/rsx-rs\"/" \
    "${DEST}/Cargo.toml" > "${tmp}"
  mv "${tmp}" "${DEST}/Cargo.toml"
fi

# Drop any accidental workspace table copied with the tree (none expected).
# Ensure no parent workspace is discovered: rsxr's Cargo.toml already has [workspace].

test -f "${DEST}/Cargo.toml"
echo "staged rsxcore ok (version field standalone)"
