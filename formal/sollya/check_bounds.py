#!/usr/bin/env python3
"""Tie a fresh Sollya run to the Rust logarithm coefficients and error constant."""

import re
import sys
from pathlib import Path


def parse_hex_float(token: str) -> float:
    token = token.strip().lower()
    sign = -1.0 if token.startswith("-") else 1.0
    if token[0] in "+-":
        token = token[1:]
    if not token.startswith("0x"):
        raise ValueError(token)
    body, exponent = token[2:].split("p")
    integer, fraction = body.split(".")
    digits = integer + fraction
    significand = int(digits, 16) / 16 ** len(fraction)
    return sign * significand * 2 ** int(exponent)


def rust_bits(source: str) -> list[int]:
    coeffs = re.findall(r"0x([0-9a-fA-F]{16}),", source)
    if len(coeffs) < 288:
        raise SystemExit(f"expected at least 288 Rust coefficients, found {len(coeffs)}")
    return [int(bits, 16) for bits in coeffs[:288]]


def rust_bound(source: str) -> float:
    match = re.search(r"LN_POLY_ABS_BOUND:\s*f64\s*=\s*([0-9.eE+-]+)", source)
    if match is None:
        raise SystemExit("LN_POLY_ABS_BOUND not found")
    return float(match.group(1))


def sollya_record(text: str) -> tuple[float, list[int]]:
    upper = None
    coeffs = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("ABS_BOUND_UPPER"):
            upper = parse_hex_float(line.split()[1])
        elif line.startswith("C "):
            token = line.split()[1]
            if token.lower().startswith("0x") and "p" not in token.lower():
                coeffs.append(int(token, 16))
    if upper is None or len(coeffs) != 288:
        raise SystemExit(f"sollya record incomplete: bound={upper} coeffs={len(coeffs)}")
    return upper, coeffs


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_bounds.py SOLLYA_OUTPUT RUST_SOURCE")
    sollya_text = Path(sys.argv[1]).read_text()
    rust_text = Path(sys.argv[2]).read_text()
    upper, sollya_coeffs = sollya_record(sollya_text)
    bound = rust_bound(rust_text)
    code_coeffs = rust_bits(rust_text)
    if sollya_coeffs != code_coeffs:
        for index, (left, right) in enumerate(zip(sollya_coeffs, code_coeffs)):
            if left != right:
                raise SystemExit(f"coefficient {index} differs: sollya {left:#x} rust {right:#x}")
        raise SystemExit("coefficient lists differ in length")
    if not (upper <= bound < 1e-15):
        raise SystemExit(f"bound {bound} does not sit above the Sollya endpoint {upper}")
    if bound > upper * 4:
        raise SystemExit(f"bound {bound} is loose relative to {upper}")
    print(f"ok bound={bound:.3e} sollya_upper={upper:.3e} coeffs=288")


if __name__ == "__main__":
    main()
