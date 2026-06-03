# Security Policy

rsx-rs is a command-line scientific tool for genomics data. It processes user-supplied FASTQ and marker tables.

## Supported versions

Only the latest release on the main branch is supported.

## Reporting a vulnerability

Please email the maintainer (rgoswami@ieee.org) with details. Do not open a public issue for security matters.

We will acknowledge within 72 hours and aim for a fix + coordinated disclosure for anything that could lead to arbitrary code execution, data corruption on trusted inputs, or supply-chain issues in the build/release artifacts.

Because the project is GPL and aimed at research use, we treat reproducible builds and provenance (pixi.lock, cargo lock, release manifests) as part of the security surface.
