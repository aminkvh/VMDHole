#!/usr/bin/env python3
"""Heteroatoms in the lining: the branch no reference run reaches.

MOLE's shipped test PDBs (1BL8, 2ACE, 1MXT) carry ZERO HETATM records - they are
stripped - and 2OAR has three gold ions. 1tqn is the only reference structure
with a real heteroatom group, its heme, and at test.xml's origin the three
tunnels never touch it: excluding HEM from the free set leaves every one of the
631 profile rows unchanged. So `IsHetAtom`'s "not a standard residue" half, and
with it the `GetResidueProperties -> null` path in both property blocks, has no
ground-truth coverage at all.

This puts a tunnel at the heme iron, where HEM must line it, and checks the
arithmetic MOLE specifies for that case: a residue absent from the property
table is skipped by the side-chain loop BEFORE the denominator is incremented,
so it appears in the lining and contributes nothing to any property. Getting
that wrong - counting it, or dividing by it - is invisible everywhere else.

Not a MOLE comparison; MOLE cannot be re-run here. It is an independent
re-derivation of the property block from the table, which is why it checks every
layer rather than only the ones holding HEM.

Usage: mole_het_check.py ENGINE_OUT.txt
"""
import sys

# TunnelPhysicoChemicalPropertyTable: hydropathy, hydrophobicity, polarity,
# logp, logd, logs, charge, ionizable, mutability. Transcribed from the C#
# independently of mole_lining.c - a shared typo is the thing this catches.
T = {
 'ALA': (1.8, 0.02, 0.0, 1.08, 1.08, 0.59, 0, 0, 100),
 'ARG': (-4.5, -0.42, 52.0, -0.08, -2.49, 1.63, 1, 1, 83),
 'ASN': (-3.5, -0.77, 3.38, -1.03, -1.03, 0.54, 0, 0, 104),
 'ASP': (-3.5, -1.04, 49.7, -0.22, -3.0, 2.63, -1, 1, 86),
 'CYS': (2.5, 0.77, 1.48, 0.84, 0.84, 0.16, 0, 0, 44),
 'GLU': (-3.5, -1.14, 49.9, 0.48, -2.12, 2.23, -1, 1, 77),
 'GLN': (-3.5, -1.1, 3.53, -0.33, -0.33, 0.13, 0, 0, 84),
 'GLY': (-0.4, -0.8, 0.0, 0.0, 0.0, 0.0, 0, 0, 50),
 'HIS': (-3.2, 0.26, 51.6, -0.01, -0.11, -0.2, 0, 0, 91),
 'ILE': (4.5, 1.81, 0.13, 2.24, 2.24, -1.85, 0, 0, 103),
 'LEU': (3.8, 1.14, 0.13, 2.08, 2.08, -1.79, 0, 0, 54),
 'LYS': (-3.9, -0.41, 49.5, 0.7, -1.91, 1.46, 1, 1, 72),
 'MET': (1.9, 1.0, 1.43, 1.48, 1.48, -0.72, 0, 0, 93),
 'PHE': (2.8, 1.35, 0.35, 2.49, 2.49, -1.81, 0, 0, 51),
 'PRO': (-1.6, -0.09, 1.58, 1.8, 1.8, -1.3, 0, 0, 58),
 'SER': (-0.8, -0.97, 1.67, -0.52, -0.52, 1.11, 0, 0, 117),
 'THR': (-0.7, -0.77, 1.66, -0.16, -0.16, 0.77, 0, 0, 107),
 'TRP': (-0.9, 1.71, 2.1, 2.59, 2.59, -2.48, 0, 0, 25),
 'TYR': (-1.3, 1.11, 1.61, 2.18, 2.18, -1.44, 0, 0, 50),
 'VAL': (4.2, 1.13, 0.13, 1.8, 1.8, -1.3, 0, 0, 98),
}
# The backbone half adds ASN's polarity, GLY's hydrophobicity/hydropathy and
# BACKBONE's logs - three different rows, deliberately.
BB = (-0.40, -0.80, 3.38, -0.86, -0.86, 0.81)
TOL = 0.00005 + 1e-9      # the engine prints properties to four decimals


def main():
    fails, layers, het_layers, het_names = 0, 0, 0, set()
    for line in open(sys.argv[1]):
        f = line.split()
        if not f or f[0] != 'L':
            continue
        layers += 1
        res = [t.split(':') for t in f[21:]]
        if any(r[0] not in T for r in res):
            het_layers += 1
            het_names |= {r[0] for r in res if r[0] not in T}

        # Re-derive the six averaged properties and the five counted ones.
        n = 0
        acc = [0.0] * 6
        charge = ion = npos = nneg = 0
        mut = 0.0
        nmut = 0
        for rn, _sq, _ch, bb, _fi in res:
            if bb == '1':
                n += 1                       # counted with NO table lookup
                acc = [a + b for a, b in zip(acc, BB)]
                continue
            if rn not in T:
                continue                     # skipped BEFORE the denominator
            v = T[rn]
            n += 1
            acc = [a + b for a, b in zip(acc, v[:6])]
            nmut += 1
            charge += v[6]
            ion += v[7]
            mut += v[8]
            npos += v[6] > 0
            nneg += v[6] < 0
        want = [a / n for a in acc] if n else [0.0] * 6
        if nmut:
            mut /= nmut
        got = [float(x) for x in (f[13], f[14], f[15], f[16], f[17], f[18])]
        for w, g, nm in zip(want, got, ('hydropathy', 'hydrophobicity',
                                        'polarity', 'logp', 'logd', 'logs')):
            if abs(w - g) > TOL:
                print('  %-52s %s' % ('layer %s/%s %s' % (f[1], f[2], nm),
                                      'FAIL (%.4f vs %.4f)' % (w, g)))
                fails += 1
        for w, g, nm in ((charge, f[9], 'charge'), (ion, f[10], 'ionizable'),
                         (npos, f[11], 'npos'), (nneg, f[12], 'nneg'),
                         (int(mut), f[19], 'mutability')):
            if int(w) != int(g):
                print('  %-52s %s' % ('layer %s/%s %s' % (f[1], f[2], nm),
                                      'FAIL (%s vs %s)' % (w, g)))
                fails += 1

    ok = het_layers > 0
    print('  %-52s %s' % ('heteroatoms reach the lining (%s)' % (
        ','.join(sorted(het_names)) or 'none'),
        'PASS (%d of %d layers)' % (het_layers, layers) if ok
        else 'FAIL - this run does not exercise the het path'))
    print('  %-52s %s' % ('properties re-derived independently, %d layers' % layers,
                          'PASS' if not fails else 'FAIL'))
    return 1 if (fails or not ok) else 0


sys.exit(main())
