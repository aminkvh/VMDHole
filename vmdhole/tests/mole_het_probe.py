#!/usr/bin/env python3
"""Why FindHetResidues' profile-sphere half never fires.

Deleting that half changes no answer on any reference, which looks exactly like
a broken transcription. It is not: MOLE's Radius is the distance from the
centreline to the nearest atom SURFACE, and InvariantKdAtomTree is keyed on atom
CENTRES, so the nearest centre sits at Radius + vdW. A search of 1.2 * Radius
reaches it only when Radius >= 5 * vdW - about 8 A, which is not a tunnel.

This asserts that geometry rather than restating it, from MOLE_HET_PROBE's
output: every probe point must have 1.2R below the nearest atom centre, and the
gap must be the nearest atom's vdW radius (1.0 - 2.31 A across MOLE's tables).

If this ever goes red, the sphere half has started to matter and the "it is
inert" reasoning in mole_lining.c no longer holds.

Usage: mole_het_probe.py PROBE.log
"""
import re
import sys

ROW = re.compile(r"HETPROBE t=\S+ R=(\S+) 1\.2R=(\S+) nearestAtomCentre=(\S+)")


def main():
    rows = [tuple(map(float, m.groups()))
            for m in map(ROW.search, open(sys.argv[1])) if m]
    bad = 0
    if not rows:
        print("  %-52s FAIL (no probe rows - was MOLE_HET_PROBE set?)" % "probe produced data")
        return 1
    reach = [r for r in rows if r[2] <= r[1]]
    gaps = [r[2] - r[0] for r in rows]
    for label, ok, detail in (
        ("%d probe points, none reached by 1.2 x Radius" % len(rows),
         not reach, "(%d reached)" % len(reach)),
        ("the gap is a vdW radius (%.3f - %.3f A)" % (min(gaps), max(gaps)),
         1.0 - 1e-9 <= min(gaps) and max(gaps) <= 2.31 + 1e-9,
         "(%.3f .. %.3f)" % (min(gaps), max(gaps))),
    ):
        print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
        if not ok:
            bad = 1
    return bad


sys.exit(main())
