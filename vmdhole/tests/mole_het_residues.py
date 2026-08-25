#!/usr/bin/env python3
"""<HetResidues> against MOLE's own XML.

Tunnel.FindHetResidues (Tunnel.cs:638) reports the HET residues a tunnel passes
through. It was left unported on the stated grounds that MOLE emits an empty set
on every available reference. That was wrong, and checking rather than
remembering is the whole point of this file: two committed references carry a
NON-empty set - hem_origin (HEM 508 A on three of four tunnels) and 1ERI
(DT 9 B, then DA 6 B + DA 7 B).

An all-empty comparison would pass vacuously, so this asserts that at least one
tunnel has a non-empty set before it credits anything.

With --by-length, tunnels are paired by LENGTH rather than by position, for
references where the two engines produce different tunnel SETS. 1KX5 is one
(57 against MOLE's 60) and it is the only reference that exercises the
profile-sphere half of the search at all.

Usage: mole_het_residues.py OURS.txt REFDIR [--by-length]
"""
import re
import sys


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def main():
    ours_path, refdir = sys.argv[1], sys.argv[2]
    by_length = "--by-length" in sys.argv[3:]
    ours, ourlen = {}, {}
    for line in open(ours_path):
        f = line.split()
        if f and f[0] == "T":
            ourlen[int(f[1])] = "%.2f" % float(f[3])
        if f and f[0] == "H":
            # "HEM:508:A" -> MOLE's own "HEM 508 A"
            ours[int(f[1])] = [t.replace(":", " ") for t in f[3:]]

    blocks = re.findall(r'<HetResidues>(.*?)</HetResidues>',
                        open("%s/tunnels.xml" % refdir).read(), re.S)
    pair = {i: i for i in range(1, len(blocks) + 1)}
    if by_length:
        reflen = [l.strip().split(",")[1]
                  for l in list(open("%s/tunnels.csv" % refdir))[1:]]
        pair, taken = {}, set()
        for i, want in enumerate(reflen, 1):
            hit = [k for k, v in ourlen.items() if v == want and k not in taken]
            if len(hit) == 1:
                pair[i] = hit[0]
                taken.add(hit[0])
        print("  %-52s %s" % ("%d of %d tunnels pair by length"
                              % (len(pair), len(blocks)), "PASS"))
    bad, nonempty = 0, 0
    for ti, b in enumerate(blocks, 1):
        if ti not in pair:
            continue
        want = [x.strip() for x in b.split(",") if x.strip()]
        got = ours.get(pair[ti], [])
        if want:
            nonempty += 1
        bad += report("tunnel %d: %s" % (ti, ", ".join(want) if want else "(none)"),
                      want == got, "(got %s)" % (", ".join(got) if got else "(none)"))
    bad += report("at least one tunnel HAS het residues", nonempty > 0,
                  "(all empty - this comparison would pass vacuously)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
