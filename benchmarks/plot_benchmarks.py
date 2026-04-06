#!/usr/bin/env python3
"""Plot benchmark results comparing C++ radsex vs Rust rsx-rs."""

import csv
from pathlib import Path
from collections import defaultdict

def load_csv(path):
    data = defaultdict(lambda: defaultdict(dict))
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            scale = row["scale"]
            cmd = row["command"]
            impl_ = row["impl"]
            t = float(row["time_secs"])
            data[scale][cmd][impl_] = t
    return data

def print_table(data):
    scales = ["small", "medium", "large"]
    commands = ["process", "freq", "depth", "distrib", "signif", "subset", "map"]

    # Header
    print(f"{'Command':<10} | {'Scale':<8} | {'C++ (s)':>8} | {'Rust (s)':>8} | {'Speedup':>8} | {'Winner':>6}")
    print("-" * 65)

    for cmd in commands:
        for scale in scales:
            if cmd not in data.get(scale, {}):
                continue
            d = data[scale][cmd]
            cpp = d.get("cpp", None)
            rust = d.get("rust", None)
            if cpp is None or rust is None:
                continue
            if rust > 0:
                speedup = cpp / rust
            else:
                speedup = float("inf")
            winner = "Rust" if speedup > 1.0 else "C++"
            marker = ">>>" if speedup > 1.5 else ("<<<" if speedup < 0.8 else "   ")
            print(f"{cmd:<10} | {scale:<8} | {cpp:>8.3f} | {rust:>8.3f} | {speedup:>7.2f}x | {winner:>5} {marker}")

    print()

    # Summary
    total_cpp = 0
    total_rust = 0
    n = 0
    for scale in scales:
        for cmd in commands:
            d = data.get(scale, {}).get(cmd, {})
            if "cpp" in d and "rust" in d:
                total_cpp += d["cpp"]
                total_rust += d["rust"]
                n += 1

    print(f"Total across {n} benchmarks: C++ {total_cpp:.3f}s, Rust {total_rust:.3f}s")
    print(f"Overall: Rust is {total_cpp/total_rust:.2f}x {'faster' if total_rust < total_cpp else 'slower'}")

    # Per-command averages
    print()
    print("Per-command average speedup (Rust/C++):")
    for cmd in commands:
        ratios = []
        for scale in scales:
            d = data.get(scale, {}).get(cmd, {})
            if "cpp" in d and "rust" in d and d["rust"] > 0:
                ratios.append(d["cpp"] / d["rust"])
        if ratios:
            avg = sum(ratios) / len(ratios)
            print(f"  {cmd:<10}: {avg:.2f}x {'faster' if avg > 1 else 'slower'}")

    # Per-scale averages
    print()
    print("Per-scale average speedup:")
    for scale in scales:
        ratios = []
        for cmd in commands:
            d = data.get(scale, {}).get(cmd, {})
            if "cpp" in d and "rust" in d and d["rust"] > 0:
                ratios.append(d["cpp"] / d["rust"])
        if ratios:
            avg = sum(ratios) / len(ratios)
            print(f"  {scale:<10}: {avg:.2f}x {'faster' if avg > 1 else 'slower'}")


if __name__ == "__main__":
    csv_path = Path("benchmarks/results/benchmark_results.csv")
    data = load_csv(csv_path)
    print_table(data)
