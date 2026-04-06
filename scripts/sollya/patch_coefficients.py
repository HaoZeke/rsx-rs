#!/usr/bin/env python3
"""
Parse Sollya erfc_coefficients.txt and patch rsx-core/src/stats.rs
with the exact f64 bit patterns.

Usage: python patch_coefficients.py [--sollya-output FILE] [--rust-source FILE]
       Or: pixi run patch-coefficients
"""

import argparse
import re
import struct
from pathlib import Path


def hex_float_to_bits(s: str) -> int:
    """Convert a hex float string (e.g. '0x1.20dd750429b6dp0') to f64 bit pattern."""
    return struct.unpack("<Q", struct.pack("<d", float.fromhex(s.strip())))[0]


def parse_sollya_output(text: str) -> dict[str, list[str]]:
    """Extract coefficient arrays from Sollya output."""
    arrays = {}
    current_name = None
    current_coeffs = []

    for line in text.splitlines():
        # Match array declaration: const ERFC_R1: [f64; 21] = [
        m = re.match(r"const\s+(\w+):\s*\[f64;\s*\d+\]\s*=\s*\[", line)
        if m:
            current_name = m.group(1)
            current_coeffs = []
            continue

        if current_name and line.strip() == "];":
            arrays[current_name] = current_coeffs
            current_name = None
            continue

        if current_name:
            # Extract hex float from line like "     0x1.20dd750429b6dp0 ,"
            m = re.search(r"(-?0x[0-9a-fA-F.]+p[+-]?\d+)", line.strip())
            if m:
                current_coeffs.append(m.group(1))

    return arrays


def generate_rust_array(name: str, coeffs: list[str]) -> str:
    """Generate a Rust const array with f64::from_bits."""
    n = len(coeffs)
    lines = [f"#[rustfmt::skip]", f"const {name}: [f64; {n}] = ["]
    for i in range(0, len(coeffs), 2):
        pair = []
        for j in range(2):
            if i + j < len(coeffs):
                bits = hex_float_to_bits(coeffs[i + j])
                pair.append(f"f64::from_bits(0x{bits:016X})")
        lines.append("    " + ", ".join(pair) + ",")
    lines.append("];")
    return "\n".join(lines)


def patch_rust_source(source: str, arrays: dict[str, list[str]]) -> str:
    """Replace coefficient arrays in Rust source with new values."""
    for name, coeffs in arrays.items():
        new_array = generate_rust_array(name, coeffs)
        # Find and replace the existing array
        pattern = rf"#\[rustfmt::skip\]\s*\nconst {name}: \[f64; \d+\] = \[.*?\];"
        source = re.sub(pattern, new_array, source, flags=re.DOTALL)
    return source


def main():
    parser = argparse.ArgumentParser(description="Patch Sollya coefficients into Rust")
    parser.add_argument(
        "--sollya-output",
        default="sollya/erfc_coefficients.txt",
        help="Sollya output file",
    )
    parser.add_argument(
        "--rust-source",
        default="../rsx-core/src/stats.rs",
        help="Rust source file to patch",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print but don't write")
    args = parser.parse_args()

    sollya_text = Path(args.sollya_output).read_text()
    arrays = parse_sollya_output(sollya_text)

    print(f"Parsed {len(arrays)} coefficient arrays:")
    for name, coeffs in arrays.items():
        print(f"  {name}: {len(coeffs)} coefficients")

    rust_source = Path(args.rust_source).read_text()
    patched = patch_rust_source(rust_source, arrays)

    if args.dry_run:
        print("\n--- Patched source (dry run) ---")
        for name in arrays:
            # Find the patched array
            m = re.search(rf"const {name}:.*?;", patched, re.DOTALL)
            if m:
                print(m.group(0))
                print()
    else:
        Path(args.rust_source).write_text(patched)
        print(f"\nPatched {args.rust_source}")


if __name__ == "__main__":
    main()
