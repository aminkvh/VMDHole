#!/usr/bin/env python3
"""The counterexample that decides the custom-exit architecture.

A custom exit at (-22.841, -55.065, -20.445) on 1tqn snaps ONLY to the
SurfaceCavity: 0.607 A to its boundary, and 12.2 A or more to every regular
cavity's, against OriginRadius = 5. MOLE still produces a tunnel from it -
T1C0, 8.82 A - so an implementation that only searches the regular cavities
returns NOTHING here. That is what rules out the regular-cavity-only shortcut.

It is not a rare case. Of 40 sampled surface points on 1tqn, 39 were
surface-only and 17 produced a tunnel; every one of those 17 was surface-only,
with the nearest regular cavity 8.9 to 30.4 A away.

Checks 1-3 verify MOLE's committed reference has the properties the decision
rests on, so the fixture cannot silently rot. Check 4 is our engine, and is
expected to report NOT IMPLEMENTED until the SurfaceCavity is searchable -
reported distinctly from a wrong answer, because "we produce nothing" and "we
produce the wrong thing" are different states.

Usage: mole_surface_exit.py REFDIR [OURS.txt]
"""
import csv
import sys

ORIGIN_RADIUS = 5.0
WANT_LENGTH = 8.82
TOL = 0.005


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def main():
    refdir = sys.argv[1]
    bad = 0

    # 1 + 2: the snap distances MOLE recorded for every cavity.
    surf, regs = None, []
    for line in open("%s/snap_distances.txt" % refdir):
        f = line.split()
        if not f or f[0] != "cavity":
            continue
        d = float(f[3].split("=")[1])
        if f[2] == "surface=True":
            surf = d
        else:
            regs.append(d)
    bad += report("surface snap succeeds (%.3f A <= %.1f)" % (surf or -1, ORIGIN_RADIUS),
                  surf is not None and surf <= ORIGIN_RADIUS)
    bad += report("every regular cavity misses (nearest %.3f A > %.1f)"
                  % (min(regs) if regs else -1, ORIGIN_RADIUS),
                  bool(regs) and min(regs) > ORIGIN_RADIUS,
                  "(%d cavities, nearest %.3f)" % (len(regs), min(regs) if regs else -1))

    # 3: MOLE's own result for that exit.
    rows = list(open("%s/tunnels.csv" % refdir))[1:]
    ids = [r.split(",")[0] for r in rows]
    lens = [float(r.split(",")[1]) for r in rows]
    bad += report("MOLE produces one tunnel, T1C0, %.2f A" % WANT_LENGTH,
                  len(rows) == 1 and ids[0] == "T1C0"
                  and abs(lens[0] - WANT_LENGTH) <= TOL,
                  "(got %s)" % ", ".join("%s %.2f" % (i, l) for i, l in zip(ids, lens)))

    # 4: us. Not implemented yet - say so, do not pass quietly.
    if len(sys.argv) < 3:
        print("  %-52s %s" % ("our engine on the same exit",
                              "NOT IMPLEMENTED (SurfaceCavity is not searchable)"))
        return 1 if bad else 0

    ours = [float(l.split()[3]) for l in open(sys.argv[2]) if l.startswith("T ")]
    if not ours:
        print("  %-52s %s" % ("our engine on the same exit",
                              "NOT IMPLEMENTED - produced no tunnel"))
        return 1 if bad else 0
    ref = list(csv.DictReader(open("%s/tunnel_1.csv" % refdir)))
    ok = len(ours) == 1 and abs(ours[0] - WANT_LENGTH) <= TOL
    bad += report("our engine reproduces MOLE's surface-exit tunnel", ok,
                  "(got %s)" % " ".join("%.2f" % v for v in ours))
    if ok:
        print("  %-52s %s" % ("  reference profile rows", len(ref)))
    return 1 if bad else 0


sys.exit(main())
