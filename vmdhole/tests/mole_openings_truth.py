#!/usr/bin/env python3
"""Cavity openings against MOLE's own, tetrahedron index for tetrahedron index.

Only checkable because the triangulation reproduces MOLE's cell order, so the
two sides mean the same thing by "tetrahedron 53844".

What it guards: MOLE takes a cavity's boundary tetrahedra as
`graph.Vertices.Where(v => v.IsBoundary)` (Cavity.cs:235), i.e. in
cavityGraph.Vertices order, not by index. That order decides MaxBy's tie
(MaxBy returns ALL maxima and Create takes [0]), Cover's recursion, and the
connected-component ids. Ordering by index instead cost two openings on 1KX5 -
cavity 1 gave 21 against 22 and cavity 3 gave 5 against 6 - and silently moved
the pivot in two more cavities without changing their count, which is the case a
count-only check would have missed.

Usage: mole_openings_truth.py OURS_OPENING_DEBUG.log MOLE_OPENINGS.txt
"""
import sys


def main():
    mole = {}
    for line in open(sys.argv[2]):
        f = line.split()
        mole[int(f[1])] = set(map(int, f[3:]))
    ours = {}
    for line in open(sys.argv[1]):
        if not line.startswith("OPENINGS"):
            continue
        f = line.split()
        ours[int(f[2])] = set(map(int, f[f.index(":") + 1:]))

    bad = 0
    for c in sorted(mole):
        m, o = mole[c], ours.get(c, set())
        if m != o:
            bad += 1
            print("  cavity %-3d MOLE %d ours %d   MOLE-only %s ours-only %s"
                  % (c, len(m), len(o), sorted(m - o), sorted(o - m)))
    tm = sum(len(v) for v in mole.values())
    to = sum(len(v) for v in ours.values())
    print("  %-52s %s" % ("%d cavities, %d openings, index for index"
                          % (len(mole), tm),
                          "PASS" if bad == 0 and tm == to else "FAIL"))
    return 1 if bad else 0


sys.exit(main())
