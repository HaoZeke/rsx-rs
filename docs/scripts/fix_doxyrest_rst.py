#!/usr/bin/env python3
"""Fix Doxyrest RST that Sphinx mishandles; rewrite thin api/index landing page."""
from __future__ import annotations

import re
import sys
from pathlib import Path

API = Path(__file__).resolve().parents[1] / "source" / "api"

INDEX_RST = """API Reference
=============

Generated C FFI documentation from ``rsx.h`` (Doxygen → Doxyrest). Prefer the
narrative :doc:`/reference/c-api` page for the contract and language notes; use
this section for symbol-level detail.

.. toctree::
   :maxdepth: 2

   global
   enum_rsx_status_t

* :doc:`global` — typedefs, enums, and entry points
* :doc:`enum_rsx_status_t` — status codes returned by fallible calls
* :ref:`genindex` — full index
"""


def fix_text(t: str) -> str:
    lines = t.splitlines(keepends=True)
    out = []
    in_toc = False
    for line in lines:
        # normalize tabs from doxyrest
        line = line.replace("\t", "   ")
        if re.match(r"^\.\. toctree::\s*$", line):
            in_toc = True
            out.append(line)
            continue
        if in_toc:
            stripped = line.lstrip()
            if stripped.startswith(":") and not line.startswith("   :"):
                out.append("   " + stripped)
                continue
            if line.strip() == "":
                out.append(line)
                continue
            # toctree entries must be indented
            if stripped and not line.startswith(" ") and not stripped.startswith("."):
                # entry like global.rst
                name = stripped.strip()
                if name.endswith(".rst"):
                    name = name[: -4]
                out.append(f"   {name}\n")
                continue
            in_toc = False
        m = re.match(r"^\|\s*:doc:`([^`]+)`\s*$", line.strip())
        if m:
            out.append(f"* :doc:`{m.group(1)}`\n")
            continue
        m = re.match(r"^\|\s*:ref:`([^`]+)`\s*$", line.strip())
        if m:
            out.append(f"* :ref:`{m.group(1)}`\n")
            continue
        out.append(line)
    t = "".join(out)
    t = re.sub(
        r"(\.\. ref-code-block::[^\n]+)\n:class:",
        r"\1\n   :class:",
        t,
    )
    t = re.sub(
        r"(\.\. code-block::[^\n]+)\n:class:",
        r"\1\n   :class:",
        t,
    )
    return t


def main() -> int:
    if not API.is_dir():
        print(f"fix_doxyrest_rst: no {API}", file=sys.stderr)
        return 0
    # Always write a usable landing page (doxyrest index is nearly empty).
    idx = API / "index.rst"
    idx.write_text(INDEX_RST, encoding="utf-8")
    print("rewrote api/index.rst")
    n = 1
    for path in sorted(API.rglob("*.rst")):
        if path.name == "index.rst":
            continue
        orig = path.read_text(encoding="utf-8")
        new = fix_text(orig)
        if new != orig:
            path.write_text(new, encoding="utf-8")
            n += 1
            print(f"fixed {path.relative_to(API.parent)}")
    print(f"fix_doxyrest_rst: {n} files touched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
