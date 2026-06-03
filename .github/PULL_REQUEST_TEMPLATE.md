## Summary
One sentence what this does and why.

## Checklist
- [ ] `cargo fmt --check && cargo clippy -D warnings` clean
- [ ] New or changed public surface has rustdoc / CLI help text
- [ ] Tests added or updated (unit / integration / precision)
- [ ] Docs updated in `docs/orgmode/` (commands, quickstart, etc.) if user-visible
- [ ] CHANGELOG.md entry added under [Unreleased] if noteworthy
- [ ] `rsx --help` and relevant subcommand help look good

## Performance / reproducibility impact
If this touches hot paths or benchmarked commands, attach before/after timings (synthetic or literature) or note that the golden outputs are unchanged.

## Related
- Fixes #...
- Paper / repro impact: ...
