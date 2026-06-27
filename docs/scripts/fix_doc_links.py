#!/usr/bin/env python3
"""Post-process exported RST: links, escaped :doc: roles, org/markdown bold leaks."""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "source"


def fix_text(t: str) -> str:
    def repl_link(m: re.Match[str]) -> str:
        target = m.group(2).strip()
        for ext in (".rst", ".org"):
            if target.endswith(ext):
                target = target[: -len(ext)]
                break
        return f":doc:`{target}`"

    t = re.sub(r"`([^\`<>]+)\s+<([^>]+?\.(?:rst|org))>`_", repl_link, t)
    # ox-rst escapes roles as :doc:\`...\` — restore
    t = re.sub(r":(doc|ref|mod|class|func|meth|attr|exc|data|const|envvar|token|option|term|eq|abbr|menuselection|guilabel|kbd|command|program|makevar|dfn|file|samp|pep|rfc|mailheader|mimetype|newsgroup|code):\\`([^`]+)\\`", r":\1:`\2`", t)
    # double-escaped angle brackets in roles
    t = t.replace("\\<", "<").replace("\\>", ">")
    # ****Label**** or **Label** left as strong+stars from org ** inside paragraphs
    t = re.sub(r"\*\*\*\*([^*]+)\*\*\*\*", r"**\1**", t)
    # Author line noise from ox-rst (duplicate of HTML theme)
    t = re.sub(r"^:Author:.*\n", "", t, flags=re.M)
    t = re.sub(r"^\.\. sectionauthor::.*\n", "", t, flags=re.M)
    t = re.sub(r"(\S)\./(\s|$)", r"\1.\2", t)
    return t


def main() -> int:
    if not ROOT.is_dir():
        print(f"fix_doc_links: missing {ROOT}", file=sys.stderr)
        return 1
    n = 0
    for path in sorted(ROOT.rglob("*.rst")):
        if "crates/" in path.as_posix() or "/xml/" in path.as_posix():
            continue
        if path.name.startswith("group__") or path.name.startswith("struct_"):
            continue
        # still fix api/ for author/escape issues
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
