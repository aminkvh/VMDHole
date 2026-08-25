#!/usr/bin/env python3
"""Compare the engine's lining layers against MOLE 2's OWN output.

The companion of mole_ground_truth.py, and it discriminates where that one
cannot: tunnels.csv carries ten heavily-aggregated scalars per tunnel, so a
misclassified backbone atom or a shifted merge boundary can cancel out inside a
length-weighted mean and still land within half a printed digit. This checks the
layers themselves - how many, where they start and end, which residues line each
and which of those are touched only at backbone, the minimum widths over each
run, the local-minimum flags, and every per-layer property.

Reference: MOLE 2's own tunnels.xml for its shipped 1tqn test, which prints the
widths to FIVE decimals (the CSVs give two) and names the lining residues.

With --by-length, tunnels are paired by LENGTH instead of by position, and
MOLE tunnels with no counterpart are reported and skipped. That exists for 1MXT,
where MOLE's own removal defect gives the two engines different tunnel SETS
(ours 8.44/11.93/25.31/3.47, MOLE's 11.93/12.04/25.31/3.47) - the three shared
ones can still be compared, and 1MXT is the only fixture carrying hydrogens, so
without this MOLE's `backboneNames` entry for "H" is never checked against MOLE.

Usage: mole_lining_truth.py OURS.txt REFDIR [--by-length]
"""
import re
import sys

GEOM_TOL = 0.000005 + 1e-12    # half of MOLE's last XML digit
PROP_TOL = 0.005 + 1e-9        # properties are printed to two decimals

# Field offsets in an "L" line; see mole_main.c's header comment.
GEOM = (("StartDistance", 3), ("EndDistance", 4), ("MinRadius", 5),
        ("MinFreeRadius", 6), ("MinBRadius", 7))
INTS = (("Charge", 9), ("Ionizable", 10), ("NumPositives", 11),
        ("NumNegatives", 12), ("Mutability", 19))
REALS = (("Hydropathy", 13), ("Hydrophobicity", 14), ("Polarity", 15),
         ("LogP", 16), ("LogD", 17), ("LogS", 18))

LAYER_RE = re.compile(
    r'<Layer ([^>]*)>\s*<Residues>(.*?)</Residues>\s*'
    r'<FlowIndices>[^<]*</FlowIndices>\s*<Properties ([^/]*)/>', re.S)
ATTR_RE = re.compile(r'(\w+)="([^"]*)"')


def residue_str(token):
    """"ASN:312:A:1[:flow]" -> MOLE's own "ASN 312 A Backbone"."""
    f = token.split(":")
    return "%s %s %s%s" % (f[0], f[1], f[2], " Backbone" if f[3] == "1" else "")


def report(label, ok, detail=""):
    print("  %-52s %s" % (label, "PASS" if ok else "FAIL " + detail))
    return 0 if ok else 1


def main():
    ours_path, refdir = sys.argv[1], sys.argv[2]
    by_length = "--by-length" in sys.argv[3:]
    ours, wprops, flows, ourlen = {}, {}, {}, {}
    for line in open(ours_path):
        f = line.split()
        if not f:
            continue
        if f[0] == "T":
            ourlen[int(f[1])] = "%.2f" % float(f[3])
        elif f[0] == "L":
            ours.setdefault(int(f[1]), []).append(f)
        elif f[0] == "W":
            wprops[int(f[1])] = f
        elif f[0] == "F":
            flows[int(f[1])] = [residue_str(t) for t in f[3:]]

    blocks = re.findall(r'<Layers>(.*?)</Layers>',
                        open("%s/tunnels.xml" % refdir).read(), re.S)
    reflen = [l.strip().split(",")[1] for l in list(open("%s/tunnels.csv" % refdir))[1:]]
    bad = 0
    # MOLE index -> the index. Positional unless --by-length, in which case each
    # of MOLE's lengths claims the one tunnel of ours that matches it; anything
    # unmatched on either side is named rather than silently dropped.
    pair = {i: i for i in range(1, len(blocks) + 1)}
    if by_length:
        pair, taken = {}, set()
        for i, want in enumerate(reflen, 1):
            hit = [k for k, v in ourlen.items() if v == want and k not in taken]
            if len(hit) == 1:
                pair[i] = hit[0]
                taken.add(hit[0])
        skipped = [reflen[i - 1] for i in range(1, len(blocks) + 1) if i not in pair]
        extra = sorted(ourlen[k] for k in ourlen if k not in taken)
        bad += report("%d of %d tunnels pair by length" % (len(pair), len(blocks)),
                      bool(pair), "(none matched)")
        if skipped or extra:
            print("  %-52s %s" % ("", "MOLE-only %s, ours-only %s"
                                  % (skipped or "-", extra or "-")))
    for ti, b in enumerate(blocks, 1):
        if ti not in pair:
            continue
        ref = LAYER_RE.findall(b)
        got = ours.get(pair[ti], [])
        if len(ref) != len(got):
            bad += report("tunnel %d: %d layers" % (ti, len(ref)), False,
                          "(engine produced %d)" % len(got))
            continue
        worst, worstwhat, msgs = 0.0, "", []
        for i, (attrs, residues, props) in enumerate(ref):
            d = dict(ATTR_RE.findall(attrs))
            p = dict(ATTR_RE.findall(props))
            g = got[i]
            for key, idx in GEOM:
                diff = abs(float(d[key]) - float(g[idx]))
                if diff > worst:
                    worst, worstwhat = diff, key
            if int(d["LocalMinimum"]) != int(g[8]):
                msgs.append("L%d LocalMinimum %s vs %s" % (i + 1, d["LocalMinimum"], g[8]))
            # "GLY 177 A Backbone" - the flag is per residue per layer, and it
            # is what splits the property sum, so compare the whole string.
            # MOLE's XML prints a layer's residues ordered by FLOW index, which
            # is not the order the layer keeps internally; sort to match.
            want = [x.strip() for x in residues.split(",") if x.strip()]
            mine = [residue_str(t)
                    for t in sorted(g[21:], key=lambda t: int(t.split(":")[4]))]
            if want != mine:
                msgs.append("L%d lining %s vs %s" % (i + 1, want, mine))
            for key, idx in INTS:
                if int(p[key]) != int(g[idx]):
                    msgs.append("L%d %s %s vs %s" % (i + 1, key, p[key], g[idx]))
            for key, idx in REALS:
                if abs(float(p[key]) - float(g[idx])) > PROP_TOL:
                    msgs.append("L%d %s %s vs %s" % (i + 1, key, p[key], g[idx]))
        want = [x.strip() for x in
                re.search(r'<ResidueFlow>(.*?)</ResidueFlow>', b, re.S).group(1).split(",")
                if x.strip()]
        if want != flows.get(pair[ti], []):
            msgs.append("ResidueFlow %s vs %s" % (want, flows.get(pair[ti])))
        ok = worst <= GEOM_TOL and not msgs
        bad += report("tunnel %d: %d layers vs MOLE's own XML" % (ti, len(ref)), ok,
                      "(%s)" % ("; ".join(msgs[:3]) if msgs
                                else "%s off by %.6f" % (worstwhat, worst)))
        if ok:
            print("  %-52s %s" % ("", "exact to MOLE's 5 printed decimals" if worst == 0
                                  else "max %s diff %.6f" % (worstwhat, worst)))

    # The tunnel-level row MOLE prints is the LAYER-WEIGHTED struct, not the
    # plain one - a port that reduces the layers correctly but averages them
    # unweighted still matches every layer above and fails here.
    ref = list(open("%s/tunnels.csv" % refdir))[1:]
    cols = (("Charge", 2), ("Ionizable", 3), ("Hydropathy", 6),
            ("Hydrophobicity", 7), ("Polarity", 8), ("LogP", 9),
            ("LogD", 10), ("LogS", 11), ("Mutability", 12))
    msgs = []
    for i, line in enumerate(ref, 1):
        if i not in pair:
            continue
        f = line.strip().split(",")
        g = wprops.get(pair[i])
        if not g:
            msgs.append("no W line for tunnel %d" % i)
            continue
        vals = dict(zip(("Id", "Length", "Charge", "Ionizable", "Hydropathy",
                         "Hydrophobicity", "Polarity", "LogP", "LogD", "LogS",
                         "Mutability"), f))
        for key, idx in cols:
            if abs(float(vals[key]) - float(g[idx])) > PROP_TOL:
                msgs.append("T%d %s %s vs %s" % (i, key, vals[key], g[idx]))
    bad += report("tunnels.csv: %d rows of weighted properties"
                  % sum(1 for i in range(1, len(ref) + 1) if i in pair),
                  not msgs, "(%s)" % "; ".join(msgs[:4]))
    return 1 if bad else 0


sys.exit(main())
