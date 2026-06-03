---
name: Bug report
about: Report unexpected behavior, crashes, or wrong results
title: ''
labels: bug
assignees: ''
---

**rsx version**
`rsx --version` (or the commit / container / conda env you used)

**Platform**
- OS: [e.g. Ubuntu 22.04, macOS 14 arm64, Windows 11]
- How built/installed: [cargo, GH release binary, pixi, pip pyrsx, ...]

**Command that failed**
```bash
rsx <subcommand> ...
```

**Input size / data**
- Number of individuals / markers (or FASTQ size)
- Public accession or link to a small reproducer if possible

**Expected behavior**

**Actual behavior** (paste error, wrong counts, performance surprise, etc.)

**Additional context**
- Was this a regression vs C++ RADSex or vs an earlier rsx version?
- Any special features (`--features mpi`, parquet, etc.)?
- Logs or `RUST_LOG=debug` output if relevant.
