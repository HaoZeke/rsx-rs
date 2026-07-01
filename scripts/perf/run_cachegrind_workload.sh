#!/usr/bin/env bash
# Fixed Cachegrind workload for rsx-rs hot-path before/after comparisons.
# Usage: run_cachegrind_workload.sh <rsx-binary> <data-dir> <label> <outfile>
set -euo pipefail
RSX_BIN="${1:?rsx binary}"
DATA="${2:?data dir with markers.tsv and popmap.tsv}"
LABEL="${3:?label}"
OUT="${4:?outfile}"
MARKERS="${DATA}/markers.tsv"
POPMAP="${DATA}/popmap.tsv"
WORKDIR="${TMPDIR:-/tmp}/rsx-cg-$$"
mkdir -p "$WORKDIR"
trap 'rm -rf "$WORKDIR"' EXIT

# Workload A: signif with Fisher (stresses fisher_exact per marker + full parse)
# Workload B: depth exact median (stresses find_median + depth parse)
# Combined: run both under one cachegrind session for a single Ir total.

echo "=== cachegrind workload label=${LABEL} ===" | tee "$OUT"
echo "host=$(hostname)" | tee -a "$OUT"
echo "date=$(date -Is)" | tee -a "$OUT"
echo "rsx_bin=${RSX_BIN}" | tee -a "$OUT"
echo "markers=${MARKERS}" | tee -a "$OUT"
echo "cmd=valgrind --tool=cachegrind --cachegrind-out-file=${WORKDIR}/cg.out ${RSX_BIN} (signif fisher + depth)" | tee -a "$OUT"

# Prefer medium; fall back handled by caller
valgrind --tool=cachegrind \
  --cachegrind-out-file="${WORKDIR}/cg.out" \
  --branch-sim=no \
  "$RSX_BIN" signif \
    -t "$MARKERS" -p "$POPMAP" -o "${WORKDIR}/signif.tsv" \
    -d 5 -G M,F --test fisher \
  2>>"$OUT"

valgrind --tool=cachegrind \
  --cachegrind-out-file="${WORKDIR}/cg_depth.out" \
  --branch-sim=no \
  "$RSX_BIN" depth \
    -t "$MARKERS" -p "$POPMAP" -o "${WORKDIR}/depth.tsv" \
  2>>"$OUT"

echo "=== cg_annotate signif ===" | tee -a "$OUT"
cg_annotate "${WORKDIR}/cg.out" 2>&1 | head -80 | tee -a "$OUT"
echo "=== cg_annotate depth ===" | tee -a "$OUT"
cg_annotate "${WORKDIR}/cg_depth.out" 2>&1 | head -80 | tee -a "$OUT"

# Extract SUMMARY lines (Ir totals)
echo "=== SUMMARY_EXTRACT ===" | tee -a "$OUT"
# cachegrind summary is in the .out file footer
grep -E 'summary:|events:' "${WORKDIR}/cg.out" | tee -a "$OUT" || true
grep -E 'summary:|events:' "${WORKDIR}/cg_depth.out" | tee -a "$OUT" || true
# Also parse from stderr I refs if present
grep -E 'I\s+refs:|D\s+refs:|I1\s+misses:|LLi\s+misses:' "$OUT" | tee -a "${OUT}.metrics" || true
echo "DONE label=${LABEL}" | tee -a "$OUT"
