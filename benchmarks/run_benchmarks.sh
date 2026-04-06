#!/usr/bin/env bash
set -euo pipefail

# Benchmark harness: compare C++ radsex vs Rust rsx-rs
# Runs each command 3 times, reports median wall-clock time

CPP_BIN="${CPP_RADSEX:-../<path>"
RUST_BIN="${RUST_RADSEX:-./target/release/rsx}"
DATA_DIR="${BENCH_DATA:-./benchmarks/data}"
RESULTS_DIR="${BENCH_RESULTS:-./benchmarks/results}"
REPS=3

mkdir -p "$RESULTS_DIR"

# Check binaries exist
if [[ ! -x "$CPP_BIN" ]]; then
    echo "ERROR: C++ radsex not found at $CPP_BIN"
    exit 1
fi
if [[ ! -x "$RUST_BIN" ]]; then
    echo "ERROR: Rust radsex not found at $RUST_BIN"
    exit 1
fi

echo "C++ binary: $CPP_BIN"
echo "Rust binary: $RUST_BIN"
echo "Data: $DATA_DIR"
echo "Repetitions: $REPS"
echo ""

# Time a command, return wall-clock seconds (median of REPS runs)
bench() {
    local label="$1"
    shift
    local times=()

    for i in $(seq 1 $REPS); do
        local t=$( { time "$@" > /dev/null 2>&1; } 2>&1 | grep real | sed 's/real\t//' )
        # Convert Xm Y.ZZZs to seconds
        local mins=$(echo "$t" | sed 's/m.*//')
        local secs=$(echo "$t" | sed 's/.*m//' | sed 's/s//')
        local total=$(echo "$mins * 60 + $secs" | bc)
        times+=("$total")
    done

    # Sort and take median
    local sorted=($(printf '%s\n' "${times[@]}" | sort -n))
    local mid=$(( ${#sorted[@]} / 2 ))
    echo "${sorted[$mid]}"
}

# CSV output
CSVFILE="$RESULTS_DIR/benchmark_results.csv"
echo "scale,command,impl,time_secs" > "$CSVFILE"

for scale in small medium large; do
    SD="$DATA_DIR/$scale"
    if [[ ! -d "$SD" ]]; then
        echo "Skipping $scale (no data)"
        continue
    fi

    TABLE="$SD/markers.tsv"
    POPMAP="$SD/popmap.tsv"
    GENOME="$SD/genome.fa"

    echo "=============================="
    echo "Scale: $scale"
    echo "=============================="

    # --- freq ---
    echo -n "  freq (C++)... "
    t_cpp=$(bench "cpp_freq" "$CPP_BIN" freq -t "$TABLE" -o "$SD/cpp_freq.tsv" -d 5)
    echo "${t_cpp}s"
    echo -n "  freq (Rust)... "
    t_rust=$(bench "rust_freq" "$RUST_BIN" freq -t "$TABLE" -o "$SD/rust_freq.tsv" -d 5)
    echo "${t_rust}s"
    echo "$scale,freq,cpp,$t_cpp" >> "$CSVFILE"
    echo "$scale,freq,rust,$t_rust" >> "$CSVFILE"

    # --- distrib ---
    echo -n "  distrib (C++)... "
    t_cpp=$(bench "cpp_distrib" "$CPP_BIN" distrib -t "$TABLE" -p "$POPMAP" -o "$SD/cpp_distrib.tsv" -d 5 -G M,F)
    echo "${t_cpp}s"
    echo -n "  distrib (Rust)... "
    t_rust=$(bench "rust_distrib" "$RUST_BIN" distrib -t "$TABLE" -p "$POPMAP" -o "$SD/rust_distrib.tsv" -d 5 -G M,F)
    echo "${t_rust}s"
    echo "$scale,distrib,cpp,$t_cpp" >> "$CSVFILE"
    echo "$scale,distrib,rust,$t_rust" >> "$CSVFILE"

    # --- signif ---
    echo -n "  signif (C++)... "
    t_cpp=$(bench "cpp_signif" "$CPP_BIN" signif -t "$TABLE" -p "$POPMAP" -o "$SD/cpp_signif.tsv" -d 5 -G M,F)
    echo "${t_cpp}s"
    echo -n "  signif (Rust)... "
    t_rust=$(bench "rust_signif" "$RUST_BIN" signif -t "$TABLE" -p "$POPMAP" -o "$SD/rust_signif.tsv" -d 5 -G M,F)
    echo "${t_rust}s"
    echo "$scale,signif,cpp,$t_cpp" >> "$CSVFILE"
    echo "$scale,signif,rust,$t_rust" >> "$CSVFILE"

    # --- depth ---
    echo -n "  depth (C++)... "
    t_cpp=$(bench "cpp_depth" "$CPP_BIN" depth -t "$TABLE" -p "$POPMAP" -o "$SD/cpp_depth.tsv")
    echo "${t_cpp}s"
    echo -n "  depth (Rust)... "
    t_rust=$(bench "rust_depth" "$RUST_BIN" depth -t "$TABLE" -p "$POPMAP" -o "$SD/rust_depth.tsv")
    echo "${t_rust}s"
    echo "$scale,depth,cpp,$t_cpp" >> "$CSVFILE"
    echo "$scale,depth,rust,$t_rust" >> "$CSVFILE"

    # --- subset ---
    echo -n "  subset (C++)... "
    t_cpp=$(bench "cpp_subset" "$CPP_BIN" subset -t "$TABLE" -p "$POPMAP" -o "$SD/cpp_subset.tsv" -d 5 -G M,F -m 3 -n 0 -N 2)
    echo "${t_cpp}s"
    echo -n "  subset (Rust)... "
    t_rust=$(bench "rust_subset" "$RUST_BIN" subset -t "$TABLE" -p "$POPMAP" -o "$SD/rust_subset.tsv" -d 5 -G M,F -m 3 -n 0 -N 2)
    echo "${t_rust}s"
    echo "$scale,subset,cpp,$t_cpp" >> "$CSVFILE"
    echo "$scale,subset,rust,$t_rust" >> "$CSVFILE"

    # --- map (skip on large, takes too long for benchmark) ---
    if [[ "$scale" != "large" ]]; then
        echo -n "  map (C++)... "
        t_cpp=$(bench "cpp_map" "$CPP_BIN" map -t "$TABLE" -p "$POPMAP" -g "$GENOME" -o "$SD/cpp_map.tsv" -d 5 -G M,F -q 0)
        echo "${t_cpp}s"
        echo -n "  map (Rust)... "
        t_rust=$(bench "rust_map" "$RUST_BIN" map -t "$TABLE" -p "$POPMAP" -g "$GENOME" -o "$SD/rust_map.tsv" -d 5 -G M,F -q 0)
        echo "${t_rust}s"
        echo "$scale,map,cpp,$t_cpp" >> "$CSVFILE"
        echo "$scale,map,rust,$t_rust" >> "$CSVFILE"
    fi

    # --- process (only for small/medium where we have FASTQ files) ---
    if [[ -d "$SD/reads" ]]; then
        echo -n "  process (C++)... "
        t_cpp=$(bench "cpp_process" "$CPP_BIN" process -i "$SD/reads" -o "$SD/cpp_process.tsv" -T 4 -d 1)
        echo "${t_cpp}s"
        echo -n "  process (Rust)... "
        t_rust=$(bench "rust_process" "$RUST_BIN" process -i "$SD/reads" -o "$SD/rust_process.tsv" -T 4 -d 1)
        echo "${t_rust}s"
        echo "$scale,process,cpp,$t_cpp" >> "$CSVFILE"
        echo "$scale,process,rust,$t_rust" >> "$CSVFILE"
    fi

    echo ""
done

echo "Results saved to $CSVFILE"
echo ""
cat "$CSVFILE"
