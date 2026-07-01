#!/usr/bin/env bash
# Measure PCA Gram path + bitset + (optional) Arrow freq on a fixed table.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RSX="${1:-$ROOT/target/release/rsx}"
DATA="${2:-$ROOT/benchmarks/data/medium}"
OUT="${3:-/dev/stdout}"
{
  echo "=== tensor-path microbench ==="
  echo "host=$(hostname) date=$(date -Is)"
  echo "rsx=$RSX data=$DATA"
  # warm
  "$RSX" pca -t "$DATA/markers.tsv" -o /tmp/rsx-pca-warm -d 1 >/dev/null 2>&1 || true
  # PCA wall times (5 runs)
  echo "--- pca wall (5x) ---"
  for i in 1 2 3 4 5; do
    /usr/bin/time -f "pca_wall_sec=%e" "$RSX" pca -t "$DATA/markers.tsv" -o "/tmp/rsx-pca-b$i" -d 1 2>&1 | tail -1
  done
  echo "--- freq wall (5x) ---"
  for i in 1 2 3 4 5; do
    /usr/bin/time -f "freq_wall_sec=%e" "$RSX" freq -t "$DATA/markers.tsv" -o "/tmp/rsx-freq-b$i.tsv" -d 5 2>&1 | tail -1
  done
  echo "--- signif chisq wall (3x, bitset path) ---"
  for i in 1 2 3; do
    /usr/bin/time -f "signif_wall_sec=%e" "$RSX" signif -t "$DATA/markers.tsv" -p "$DATA/popmap.tsv" -o "/tmp/rsx-sig-b$i.tsv" -d 5 -G M,F --test chisq 2>&1 | tail -1
  done
} | tee "$OUT"
