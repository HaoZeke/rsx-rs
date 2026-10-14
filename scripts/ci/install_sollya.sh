#!/bin/sh
# Build Sollya 8.0 if it is not already on PATH. Used by the formal workflow.
set -eu
if command -v sollya >/dev/null 2>&1; then
  exit 0
fi
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential pkg-config libgmp-dev libmpfr-dev libmpfi-dev libfplll-dev libxml2-dev ca-certificates
workdir=$(mktemp -d)
curl -fsSL -o "$workdir/sollya-8.0.tar.gz" https://www.sollya.org/releases/sollya-8.0/sollya-8.0.tar.gz
tar -xzf "$workdir/sollya-8.0.tar.gz" -C "$workdir"
(
  cd "$workdir/sollya-8.0"
  ./configure --prefix=/usr/local
  make -j"$(nproc)"
  sudo make install
)
sudo ldconfig
sollya -v || true
