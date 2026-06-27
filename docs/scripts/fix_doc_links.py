#!/usr/bin/env python3
"""Rewrite accidental `label <page.rst|org>`_ links to :doc:`…` in Sphinx RST.

Run after org→RST export (mkrst). Mirrors LODE readcon-core so lychee does not
see file:// …/foo.rst links in built HTML.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "source"


def fix_text(t: str) -> str:
    def repl(m: re.Match[str]) -> str:
        target = m.group(2).strip()
        # Keep Sphinx doc path relative (drop extension); preserve ../ segments.
        stem = target
        for ext in (".rst", ".org"):
            if stem.endswith(ext):
                stem = stem[: -len(ext)]
                break
        # Absolute-ish paths under source/ → doc name from ROOT
        p = Path(stem)
        if not stem.startswith("."):
            # e.g. reference/bindings or commands
            return f":doc:`{stem}`"
        # relative ../reference/bindings → resolve later as written
        return f":doc:`{stem}`"

    t = re.sub(r"`([^\`<>]+)\s+<([^>]+?\.(?:rst|org))>`_", repl, t)
    # stray "./" after period from some ox-rst exports
    t = re.sub(r"(\S)\./(\s|$)", r"\1.\2", t)
    return t


def main() -> int:
    if not ROOT.is_dir():
        print(f"fix_doc_links: missing {ROOT}", file=sys.stderr)
        return 1
    n = 0
    for path in sorted(ROOT.rglob("*.rst")):
        # skip generated crate API dumps (often huge; no narrative links)
        if "crates/" in path.as_posix() or "/xml/" in path.as_posix():
            continue
        if path.name.startswith("group__") or path.name.startswith("struct_"):
            continue
        orig = path.read_text(encoding="utf-8")
        new = fix_text(orig)
        if new != orig:
            path.write_text(new, encoding="utf-8")
            n += 1
            print(f"fixed {path.relative_to(ROOT)}")
    print(f"fix_doc_links: {n} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
