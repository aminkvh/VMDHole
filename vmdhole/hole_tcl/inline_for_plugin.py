#!/usr/bin/env python3
"""Bundle vmdhole/hole_tcl/ into one self-contained Tcl block for the plugin.

Expands hole.tcl's `source` lines in place, so evaluation order is exactly the
standalone port's, and drops the CLI auto-run guard. The result defines the
`hole::` namespace and executes nothing at source time.

    python3 vmdhole/hole_tcl/inline_for_plugin.py > /tmp/hole_inlined.tcl

Verified: the bundle produces a .sph byte-identical to the reference-build hole (see vmdhole/hole_tcl/README.md)
(1GRM, cpoint 0 0 0, cvect 0 0 1, sample 0.25, endrad 22, seed 1).
"""
import os
import re
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    with open(os.path.join(BASE, name), encoding="utf-8") as fh:
        return fh.read()


def expand(match):
    line = match.group(0)
    found = re.search(r"(\w+)\.tcl", line)
    if not found:
        return line
    name = found.group(1) + ".tcl"
    return (f"# ---- begin {name} (inlined from vmdhole/hole_tcl/) ----\n"
            f"{load(name)}\n# ---- end {name} ----")


def bundle():
    src = re.sub(r"^source \[file join .*\]$", expand, load("hole.tcl"), flags=re.M)
    src = re.sub(r"\nif \{\[info exists argv0\].*?\n\s*hole::main \$argv\n\}\n",
                 "\n", src, flags=re.S)
    if "hole::main $argv" in src:
        raise SystemExit("CLI auto-run guard not removed")
    if "source [file join" in src:
        raise SystemExit("unexpanded source remains")
    return src


if __name__ == "__main__":
    sys.stdout.write(bundle())
