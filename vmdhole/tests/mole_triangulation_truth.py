#!/usr/bin/env python3
"""The DH triangulation against MOLE's OWN cells - order included.

Every other triangulation check in this suite compares a cell SET, and any
correct Delaunay code passes those. This one compares the order of the four
vertices WITHIN each cell, and the order of the cells themselves, against a dump
taken out of MOLE. Neither is derivable from the geometry: they are the history
of MOLE's incremental insertion, which is why mole_dh.c exists.

It matters because the first control-path sample sits on a circumcentre
equidistant from its four atoms, so which vertex is first decides which residue
sorts first, which decides a lining layer boundary. Before this, 1BL8 reported
27 30 28 29 layers where MOLE reports 26 29 27 28.

The reference keeps the first 500 cells verbatim, for diagnosis, and a sha256
over all of them, so a regression past cell 500 cannot hide.

Usage: mole_triangulation_truth.py OURS.cells REFERENCE.cells
"""
import hashlib
import sys


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def load(path):
    head, rows = {}, []
    for line in open(path):
        if line.startswith("#"):
            f = line[1:].split()
            if len(f) == 2 and f[0] in ("cells", "sha256"):
                head[f[0]] = f[1]
            continue
        line = line.strip()
        if line:
            # "-0.0000" and "0.0000" are the same centroid: C's printf keeps the
            # sign of a small negative, .NET's formatter does not. Normalised
            # here, on BOTH sides, rather than in either producer.
            rows.append(" ".join("0.0000" if f == "-0.0000" else f
                                 for f in line.split()))
    return head, rows


def main():
    ours_path, ref_path = sys.argv[1], sys.argv[2]
    _, ours = load(ours_path)
    head, ref = load(ref_path)
    bad = 0

    want_n = int(head["cells"])
    bad += report("cell count %d (MOLE's)" % want_n, len(ours) == want_n,
                  "(got %d)" % len(ours))

    n = min(len(ref), len(ours))
    first_bad = next((i for i in range(n) if ours[i] != ref[i]), None)
    bad += report("first %d cells identical, vertex order included" % n,
                  first_bad is None,
                  "(cell %s: ours %s / MOLE %s)"
                  % (first_bad, ours[first_bad] if first_bad is not None else "",
                     ref[first_bad] if first_bad is not None else ""))

    dig = hashlib.sha256("\n".join(ours).encode()).hexdigest()
    bad += report("digest over ALL cells matches MOLE's", dig == head["sha256"],
                  "(got %s, want %s)" % (dig[:16], head["sha256"][:16]))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
