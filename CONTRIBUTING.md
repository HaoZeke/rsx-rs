# Contributing to rsx-rs

Thank you for your interest! rsx-rs is a performance-critical scientific tool; we value clear, tested, reproducible changes.

## Development environment

We use [pixi](https://pixi.sh) for reproducible environments across Rust + Python + docs + benchmarks.

```bash
pixi install -e dev
pixi run -e dev check          # fmt + clippy -D + test
pixi run -e dev build
```

For Python bindings work:

```bash
pixi run -e python build-python
pixi run -e python test-python
```

See `pixi.toml` for all environments and tasks.

## Code style & gates

- `cargo fmt --check`
- `cargo clippy -- -D warnings`
- All tests must pass (`cargo test -p rsx-core`)
- New commands or public API surface need docs + a test (unit or integration in `rsx-core/tests/`).
- CLI changes should update the org docs under `docs/orgmode/` (they are the source of truth; Sphinx + rustdocgen export).

## Adding / modifying a command

1. Implement in `rsx-core/src/commands/`.
2. Wire in `rsx-cli/src/main.rs` (clap) and the Python FFI if applicable.
3. Add a golden or precision test.
4. Document the flags + semantics in `docs/orgmode/reference/commands.org`.
5. Update quickstart/examples if user-visible.
6. Run the full benchmark smoke if performance-sensitive.

## Documentation

The primary docs live in orgmode (`docs/orgmode/`). Edit there, then `pixi run -e docs mkrst && pixi run -e docs site` (or the CI does it).

## Releasing

- Update CHANGELOG (move Unreleased → versioned section).
- Bump versions via workspace.package (root Cargo.toml) + pyproject + conf.py.
- Tag `vX.Y.Z`; cargo-dist + pypi workflows produce the artifacts.
- The companion `rsx_bmc_repro` (in the reproducibility collection, snakemake + builder profile shape) + rsx-rs `repro/` orgs are the materials for the accompanying paper.

## Questions / issues

Use the issue templates. For performance or biology claims, please include `rsx --version`, platform, exact command + input size, and (if possible) a small reproducer.

## License

By contributing you agree your work is licensed under the same GPL-3.0-or-later terms as the project.
