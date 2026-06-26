# Contributing to rsxr

rsxr is the R package under `rsx-r/` in the
[rsx-rs](https://github.com/HaoZeke/rsx-rs) monorepo. Please open issues
and pull requests on that repository.

## Development

``` bash
cd rsx-r
pixi install
sh tools/stage_rsxcore.sh .   # or let configure stage ../rsxcore
pixi run document
pixi run test
pixi run check
pixi run test-pkgcheck-workflow   # structural gate for the Actions YAML
```

System requirements: Rust (`cargo`/`rustc`, see `Config/rsxr/MSRV` in
`DESCRIPTION`), and on Linux the zlib / bzip2 / xz development
libraries.

## Code of conduct

See the repository root
[CODE_OF_CONDUCT.md](https://github.com/HaoZeke/rsx-rs/blob/main/CODE_OF_CONDUCT.md).

## License

Contributions are accepted under the same GPL-3 terms as the package
(`LICENSE` / `LICENSE.md`).
