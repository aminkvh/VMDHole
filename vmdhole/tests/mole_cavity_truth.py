#!/usr/bin/env python3
"""Compare the engine's cavities against MOLE 2's own cavities.xml.

The cavity is the space a tunnel is found in, and MOLE exports it as its own
file: Volume, Depth, DepthLength, the boundary/inner residue split, and a
physico-chemical block for each. The properties here use a DIFFERENT rule from
the per-layer ones - CalculateResidueProperties' residues overload has every
residue contribute the backbone constants as well as its own side-chain values,
and its second loop increments the denominator with no table lookup at all. A
port that reuses the layer rule matches every tunnel and gets every cavity
wrong, which is exactly why this check exists separately.

Usage: mole_cavity_truth.py OURS.txt MOLE_CAVITIES.xml
"""
import re
import sys

VOL_TOL = 0.0005 + 1e-9      # MOLE prints Volume to three decimals
DEP_TOL = 1e-6               # DepthLength is G17; ours accumulates differently
PROP_TOL = 0.005 + 1e-9      # properties are printed to two decimals

CAV_RE = re.compile(r'<Cavity ([^>]*)>(.*?)</Cavity>', re.S)
ATTR_RE = re.compile(r'(\w+)="([^"]*)"')
SIDE_RE = (("Boundary", "VB"), ("Inner", "VI"))
INTS = (("Charge", 2), ("Ionizable", 3), ("NumPositives", 4), ("NumNegatives", 5),
        ("Mutability", 12))
REALS = (("Hydropathy", 6), ("Hydrophobicity", 7), ("Polarity", 8),
         ("LogP", 9), ("LogD", 10), ("LogS", 11))


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def main():
    ours_v, ours_side = {}, {}
    for line in open(sys.argv[1]):
        f = line.split()
        if not f:
            continue
        if f[0] == "V":
            ours_v[int(f[1])] = f
        elif f[0] in ("VB", "VI"):
            ours_side[(f[0], int(f[1]))] = f

    ref = CAV_RE.findall(open(sys.argv[2]).read())
    if not ref:
        return report("MOLE cavities.xml has no cavities", False)

    # MOLE exports every cavity; we emit only the ones that passed the depth
    # filters, so match on Id and check the ones we produced.
    bad, checked = 0, 0
    msgs = []
    for attrs, body in ref:
        d = dict(ATTR_RE.findall(attrs))
        # MOLE numbers CHANNELS and VOIDS separately, so Id="1" appears twice
        # in this file. Only channels can hold tunnels and only channels are
        # what the engine emits, so voids are skipped rather than collided with.
        if d["Type"] != "Cavity":
            continue
        cid = int(d["Id"])
        if cid not in ours_v:
            continue
        checked += 1
        g = ours_v[cid]
        if g[2] != d["Type"]:
            msgs.append("C%d Type %s vs %s" % (cid, d["Type"], g[2]))
        if abs(float(d["Volume"]) - float(g[3])) > VOL_TOL:
            msgs.append("C%d Volume %s vs %s" % (cid, d["Volume"], g[3]))
        if int(d["Depth"]) != int(g[4]):
            msgs.append("C%d Depth %s vs %s" % (cid, d["Depth"], g[4]))
        if abs(float(d["DepthLength"]) - float(g[5])) > DEP_TOL:
            msgs.append("C%d DepthLength %s vs %s" % (cid, d["DepthLength"], g[5]))
        for tag, key in SIDE_RE:
            m = re.search(r'<%s>\s*<Residues>(.*?)</Residues>\s*<Properties ([^/]*)/>'
                          % tag, body, re.S)
            if not m:
                msgs.append("C%d no <%s> in the reference" % (cid, tag))
                continue
            want = [x.strip() for x in m.group(1).split(",") if x.strip()]
            p = dict(ATTR_RE.findall(m.group(2)))
            gs = ours_side.get((key, cid))
            if gs is None:
                msgs.append("C%d no %s line" % (cid, key))
                continue
            mine = ["%s %s %s" % tuple(t.split(":")) for t in gs[14:]]
            if want != mine:
                extra = sorted(set(mine) - set(want))
                miss = sorted(set(want) - set(mine))
                msgs.append("C%d %s residues: %d vs %d (missing %s, extra %s)"
                            % (cid, tag, len(want), len(mine),
                               miss[:3] or "-", extra[:3] or "-"))
            for name, idx in INTS:
                if int(p[name]) != int(gs[idx]):
                    msgs.append("C%d %s %s %s vs %s" % (cid, tag, name, p[name], gs[idx]))
            for name, idx in REALS:
                if abs(float(p[name]) - float(gs[idx])) > PROP_TOL:
                    msgs.append("C%d %s %s %s vs %s" % (cid, tag, name, p[name], gs[idx]))

    bad += report("cavities match MOLE's own cavities.xml (%d checked)" % checked,
                  checked > 0 and not msgs,
                  "(%s)" % "; ".join(msgs[:4]) if msgs else "(none matched by Id)")
    return 1 if bad else 0


sys.exit(main())
