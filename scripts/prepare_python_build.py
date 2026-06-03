#!/usr/bin/env python3
"""Prepare Python extension artifacts for a maturin build."""

from pathlib import Path
from time import time_ns


def main() -> None:
    stamp = str(time_ns())
    artifact_paths = [
        Path("target/maturin/libpyrsx.so"),
        Path("target/release/libpyrsx.so"),
        Path("target/release/deps/libpyrsx.so"),
    ]
    for path in artifact_paths:
        if path.exists() and path.stat().st_size == 0:
            path.rename(path.with_name(f"{path.name}.stale-{stamp}"))


if __name__ == "__main__":
    main()
