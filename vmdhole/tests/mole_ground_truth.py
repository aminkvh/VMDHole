#!/usr/bin/env python3
"""Compare the engine's profiles against MOLE 2's OWN output, row by row.

Every other test in this suite compares us against us: the Tcl against the C,
the C against an earlier C. This one is the only check that the answer is
MOLE's. The reference is MOLE 2's own CSV export for its shipped 1tqn test
(vmdhole/tests/fixtures/mole_reference/test_xml/), produced by MOLE itself, and it
carries the columns that nothing else pins - FreeRadius and BRadius.

MOLE prints %.3f, so agreement to 0.0005 IS exact at their printed precision;
anything larger is a real difference. The engine samples at GetProfile(8), which
is the density MOLE exported at, so rows correspond one to one.

Usage: mole_ground_truth.py OURS.txt REFDIR
"""
import csv
import sys

TOL = 0.0005 + 1e-9      # half of MOLE's last printed digit
COLS = ("Radius", "FreeRadius", "BRadius", "X", "Y", "Z")


def main():
    ours_path, refdir = sys.argv[1], sys.argv[2]
    # How many tunnels the reference has. Defaults to 1tqn's 3 because that is
    # the fixture this started as; any other reference must say.
    want = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    ours = {}
    for line in open(ours_path):
        f = line.split()
        if f and f[0] == "P":
            ours.setdefault(int(f[1]), []).append([float(v) for v in f[2:8]])

    if not ours:
        print("  %-52s %s" % ("engine produced no tunnels", "FAIL"))
        return 1

    bad = 0
    for tid in sorted(ours):
        try:
            ref = list(csv.DictReader(open("%s/tunnel_%d.csv" % (refdir, tid))))
        except OSError:
            print("  %-52s %s" % ("no MOLE reference for tunnel %d" % tid, "FAIL"))
            bad += 1
            continue
        o = ours[tid]
        if len(o) != len(ref):
            print("  %-52s %s" % ("tunnel %d row count = %d" % (tid, len(ref)),
                                  "FAIL (ours %d)" % len(o)))
            bad += 1
        n = min(len(o), len(ref))
        worst, worstcol = 0.0, ""
        for i in range(n):
            got = dict(zip(("X", "Y", "Z", "Radius", "FreeRadius", "BRadius"), o[i]))
            for k in COLS:
                d = abs(got[k] - float(ref[i][k]))
                if d > worst:
                    worst, worstcol = d, k
        ok = worst <= TOL
        if not ok:
            bad += 1
        print("  %-52s %s" % (
            "tunnel %d: %d rows vs MOLE's own export" % (tid, n),
            "PASS (max %s diff %.4f)" % (worstcol, worst) if ok
            else "FAIL (%s off by %.4f)" % (worstcol, worst)))

    print("  %-52s %s" % ("MOLE reports %d tunnels here" % want,
                          "PASS" if len(ours) == want else "FAIL (we report %d)" % len(ours)))
    if len(ours) != want:
        bad += 1
    return 1 if bad else 0


sys.exit(main())
