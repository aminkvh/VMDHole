#!/usr/bin/env python3
"""Re-inline vmdhole/hole_tcl/ into vmdhole.tcl's sentinel-delimited engine region.

The plugin ships as ONE script, so the engine lives inside vmdhole.tcl rather
than beside it. Editing anything under vmdhole/hole_tcl/ therefore has no effect until
this runs:

    python3 vmdhole/hole_tcl/sync_into_plugin.py vmdhole/vmdhole.tcl

Replaces only the text between the two sentinels, leaves everything else byte
for byte, and refuses rather than guessing if the region is not exactly one
well-formed pair.
"""
import os
import sys

import inline_for_plugin

BEGIN = "# ===== BEGIN INLINED HOLE PURE-TCL ENGINE ====="
END = "# ===== END INLINED HOLE PURE-TCL ENGINE ====="


def main(path):
    with open(path, encoding="utf-8") as fh:
        lines = fh.read().split("\n")
    b = [i for i, ln in enumerate(lines) if ln.strip() == BEGIN]
    e = [i for i, ln in enumerate(lines) if ln.strip() == END]
    if len(b) != 1 or len(e) != 1 or e[0] <= b[0] + 1:
        raise SystemExit(
            f"expected exactly one well-formed sentinel pair in {path}, "
            f"found {len(b)} BEGIN and {len(e)} END")
    body = inline_for_plugin.bundle().rstrip("\n").split("\n")
    out = lines[:b[0] + 1] + body + lines[e[0]:]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(out))
    print(f"re-inlined {len(body)} lines into {path}")


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main(sys.argv[1] if len(sys.argv) > 1 else "vmdhole/vmdhole.tcl")
