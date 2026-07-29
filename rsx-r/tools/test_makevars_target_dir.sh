#!/bin/sh
# Assert that R package builds keep Cargo output beside the package sources.
set -eu

ROOT="$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)"
MAKEVARS_IN="$ROOT/rsx-r/src/Makevars.in"

grep -Fq 'RUST_TARGET_DIR = $(CURDIR)/rust/target' "$MAKEVARS_IN"
grep -Fq 'LIBDIR = $(RUST_TARGET_DIR)/release' "$MAKEVARS_IN"
grep -Fq 'export CARGO_TARGET_DIR="$(RUST_TARGET_DIR)"' "$MAKEVARS_IN"
! grep -Fq 'LIBDIR = ./rust/target/release' "$MAKEVARS_IN"

echo "test_makevars_target_dir.sh: all assertions passed"
