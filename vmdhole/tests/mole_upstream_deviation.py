#!/usr/bin/env python3
"""The ONE place our engine deliberately disagrees with MOLE, pinned.

Everywhere else the bar is exact agreement with MOLE's own output. On 1MXT it is
not, and that is a decision rather than a defect: MOLE's similarity filter picks
the 8.44 A tunnel as the survivor of its group and its removal loop then deletes
that very tunnel, leaving a 12.04 A look-alike (value-based List.Remove where
identity is required - see vmdhole/tests/fixtures/mole_reference/upstream_bug/).

We emit what MOLE's filter selected. Reviewed and accepted: report the defect
upstream, do not reproduce it, and pin the deviation here so that

  - nobody "fixes" our 8.44 A to 12.04 A believing they are improving parity, and
  - if upstream ever repairs it, this check fails and tells us to revisit.

Usage: mole_upstream_deviation.py OURS.txt
"""
import sys

# MOLE 2 e0df21c on 1MXT, default parameters, auto origins.
MOLE_SAYS = [11.93, 12.04, 25.31, 3.47]
# What MOLE's own filter selects, and what this emit: 8.44 in place of 12.04.
WE_SAY = [8.44, 11.93, 25.31, 3.47]
TOL = 0.005


def main():
    got = []
    for line in open(sys.argv[1]):
        f = line.split()
        if f and f[0] == "T":
            got.append(float(f[3]))
    got.sort()
    want = sorted(WE_SAY)

    ok = len(got) == len(want) and all(abs(a - b) <= TOL for a, b in zip(got, want))
    print("  %-52s %s" % (
        "1MXT: we emit MOLE's filter-selected tunnel",
        "PASS (%s)" % " ".join("%.2f" % v for v in got) if ok
        else "FAIL (got %s, expected %s)" % (
            " ".join("%.2f" % v for v in got), " ".join("%.2f" % v for v in want))))
    if not ok:
        # Distinguish "someone made us match MOLE" from "something else broke".
        if len(got) == len(MOLE_SAYS) and all(
                abs(a - b) <= TOL for a, b in zip(got, sorted(MOLE_SAYS))):
            print("  %-52s %s" % (
                "  -> this now MATCHES MOLE's buggy output", "was that deliberate?"))
        return 1
    print("  %-52s %s" % (
        "  (MOLE exports %s - upstream defect, not ours)"
        % " ".join("%.2f" % v for v in sorted(MOLE_SAYS)), ""))
    return 0


sys.exit(main())
