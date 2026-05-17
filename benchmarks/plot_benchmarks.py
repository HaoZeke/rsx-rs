#!/usr/bin/env python3
"""Summarize and plot benchmark results comparing C++ RADSex and Rust rsx."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


SCALES = ["small", "medium", "large"]
COMMANDS = ["process", "freq", "depth", "distrib", "signif", "subset", "map"]


def load_csv(path: Path) -> dict[str, dict[str, dict[str, float]]]:
    data: dict[str, dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        required = {"scale", "command", "impl", "time_secs"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            scale = row["scale"]
            cmd = row["command"]
            impl = row["impl"]
            data[scale][cmd][impl] = float(row["time_secs"])
    return data


def iter_pairs(data: dict[str, dict[str, dict[str, float]]]):
    for cmd in COMMANDS:
        for scale in SCALES:
            row = data.get(scale, {}).get(cmd, {})
            cpp = row.get("cpp")
            rust = row.get("rust")
            if cpp is not None and rust is not None:
                yield cmd, scale, cpp, rust


def print_table(data: dict[str, dict[str, dict[str, float]]]) -> str:
    lines = [
        f"{'Command':<10} | {'Scale':<8} | {'C++ (s)':>8} | {'Rust (s)':>8} | {'Speedup':>8} | {'Winner':>6}",
        "-" * 65,
    ]

    total_cpp = 0.0
    total_rust = 0.0
    count = 0
    ratios_by_cmd: dict[str, list[float]] = defaultdict(list)
    ratios_by_scale: dict[str, list[float]] = defaultdict(list)

    for cmd, scale, cpp, rust in iter_pairs(data):
        speedup = cpp / rust if rust > 0 else float("inf")
        winner = "Rust" if speedup > 1.0 else "C++"
        marker = ">>>" if speedup > 1.5 else ("<<<" if speedup < 0.8 else "")
        lines.append(
            f"{cmd:<10} | {scale:<8} | {cpp:>8.3f} | {rust:>8.3f} | "
            f"{speedup:>7.2f}x | {winner:>5} {marker}"
        )
        total_cpp += cpp
        total_rust += rust
        count += 1
        ratios_by_cmd[cmd].append(speedup)
        ratios_by_scale[scale].append(speedup)

    lines.append("")
    lines.append(f"Total across {count} benchmarks: C++ {total_cpp:.3f}s, Rust {total_rust:.3f}s")
    if total_rust > 0:
        direction = "faster" if total_rust < total_cpp else "slower"
        lines.append(f"Overall: Rust is {total_cpp / total_rust:.2f}x {direction}")

    lines.append("")
    lines.append("Per-command average speedup:")
    for cmd in COMMANDS:
        ratios = ratios_by_cmd.get(cmd, [])
        if ratios:
            avg = sum(ratios) / len(ratios)
            direction = "faster" if avg > 1 else "slower"
            lines.append(f"  {cmd:<10}: {avg:.2f}x {direction}")

    lines.append("")
    lines.append("Per-scale average speedup:")
    for scale in SCALES:
        ratios = ratios_by_scale.get(scale, [])
        if ratios:
            avg = sum(ratios) / len(ratios)
            direction = "faster" if avg > 1 else "slower"
            lines.append(f"  {scale:<10}: {avg:.2f}x {direction}")

    summary = "\n".join(lines)
    print(summary)
    return summary


def write_speedup_svg(data: dict[str, dict[str, dict[str, float]]], output: Path) -> None:
    rows = [(cmd, scale, cpp / rust) for cmd, scale, cpp, rust in iter_pairs(data) if rust > 0]
    width = 1120
    left = 180
    top = 28
    row_h = 24
    plot_w = 820
    height = top * 2 + row_h * len(rows) + 34
    max_speedup = max([ratio for _, _, ratio in rows] + [1.0])

    def x_pos(value: float) -> float:
        return left + (value / max_speedup) * plot_w

    palette = {"small": "#2b6cb0", "medium": "#2f855a", "large": "#b7791f"}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;font-size:13px}.axis{stroke:#444;stroke-width:1}.tick{stroke:#ddd;stroke-width:1}</style>',
        f'<line class="axis" x1="{left}" y1="{top - 8}" x2="{left}" y2="{height - 38}"/>',
    ]

    for tick in range(1, int(max_speedup) + 2):
        x = x_pos(float(tick))
        lines.append(f'<line class="tick" x1="{x:.1f}" y1="{top - 8}" x2="{x:.1f}" y2="{height - 38}"/>')
        lines.append(f'<text x="{x - 8:.1f}" y="{height - 16}">{tick}x</text>')

    for idx, (cmd, scale, ratio) in enumerate(rows):
        y = top + idx * row_h
        label = f"{cmd} ({scale})"
        bar_w = max(1.0, x_pos(ratio) - left)
        color = palette.get(scale, "#4a5568")
        lines.append(f'<text x="16" y="{y + 15}">{label}</text>')
        lines.append(f'<rect x="{left}" y="{y}" width="{bar_w:.1f}" height="16" fill="{color}"/>')
        lines.append(f'<text x="{left + bar_w + 8:.1f}" y="{y + 13}">{ratio:.2f}x</text>')

    lines.append("</svg>")
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="benchmarks/results/benchmark_results.csv", type=Path)
    parser.add_argument("--output", default="docs/figures", type=Path)
    args = parser.parse_args()

    data = load_csv(args.input)
    summary = print_table(data)

    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "benchmark_summary.txt").write_text(summary + "\n")
    write_speedup_svg(data, args.output / "benchmark_speedups.svg")
    print(f"\nWrote {args.output / 'benchmark_summary.txt'}")
    print(f"Wrote {args.output / 'benchmark_speedups.svg'}")


if __name__ == "__main__":
    main()
