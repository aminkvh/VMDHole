#!/usr/bin/env python3
"""StrictInterior in the C engine, against MOLE's own numbers.

MOLE's command line cannot set StrictInterior: the property is declared once in
ComplexParameters.cs and never assigned, never written by ToXml, never read by
FromXml. The reference values below were therefore obtained by rebuilding MOLE
with it forced true, into a separate tree so the validated oracle binary was
left alone.

The flag makes MaxClearance the maximum over the ADJACENT EDGES' Clearance, and
every one of those is still zero at that point because the line computing them
(ComplexComputation.cs:345/348, under the author's own comment "if strict, need
to update the edges 1st!") is commented out. Zero is below 2 *
InteriorThreshold, so the interior removal takes every surviving tetrahedron.

An empty result is also what a crash or a broken build produces, so the checks
that decide anything are the non-empty ones: the SurfaceCavity is snapshotted
BEFORE the interior removal, so it survives - and it GROWS, because MaxClearance
= 0 also reduces the probe peel to its second term. Its volume has to be right.

Usage: mole_strict_interior.py NORMAL.err NORMAL.txt STRICT.err STRICT.txt
                               SEXIT_NORMAL.txt SEXIT_STRICT.txt
"""
import re
import sys

# MOLE's own 1BL8 output: the Cavity and Void entries in cavities.xml, and the
# MolecularSurface entry's Volume, from the stock build and the forced rebuild.
WANT = {0: (4, 9, "24837.879"), 1: (0, 0, "26633.820")}
# MOLE's own 1tqn surface-exit run, same two builds.
WANT_EXIT = {0: 1, 1: 0}


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def counts(err_path):
    """channels, voids and the SurfaceCavity volume from one engine run."""
    text = open(err_path).read()
    m = re.search(r"(\d+) channels, (\d+) voids", text)
    v = re.search(r"volume=([0-9.]+)", text)
    if not m:
        raise SystemExit("no cavity counts in %s - engine failed?" % err_path)
    if not v:
        raise SystemExit("no SURFACE line in %s - run with MOLE_SURFACE_DEBUG=1"
                         % err_path)
    return int(m.group(1)), int(m.group(2)), v.group(1)


def tunnels(out_path):
    return sum(1 for line in open(out_path) if line.startswith("T "))


def main():
    nerr, nout, serr, sout, xnorm, xstrict = sys.argv[1:7]
    bad = 0
    for strict, err, out in ((0, nerr, nout), (1, serr, sout)):
        wch, wvd, wvol = WANT[strict]
        ch, vd, vol = counts(err)
        bad += report("StrictInterior=%d: %d channels, %d voids (MOLE's)"
                      % (strict, wch, wvd),
                      ch == wch and vd == wvd, "(got %d/%d)" % (ch, vd))
        bad += report("StrictInterior=%d: SurfaceCavity volume %s (MOLE's)"
                      % (strict, wvol), vol == wvol, "(got %s)" % vol)
        want_t = 8 if strict == 0 else 0
        got_t = tunnels(out)
        bad += report("StrictInterior=%d: %d tunnels on 1BL8 (MOLE's)"
                      % (strict, want_t), got_t == want_t, "(got %d)" % got_t)

    # The surface-exit run reaches the SurfaceCavity directly, so it is the case
    # where "strict still finds something" would be most plausible. MOLE finds
    # nothing: with no cavities there is no origin to start from.
    for strict, path in ((0, xnorm), (1, xstrict)):
        got = tunnels(path)
        bad += report("StrictInterior=%d: %d surface-exit tunnel(s) (MOLE's)"
                      % (strict, WANT_EXIT[strict]),
                      got == WANT_EXIT[strict], "(got %d)" % got)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
