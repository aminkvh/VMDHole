#!/usr/bin/env python3
"""Emit a PDB's ATOM/HETATM records in FILE ORDER, for the MOLE port's oracle test.

The sibling of extract_atoms.py, which does the same for mmCIF. MOLE's own test
suite ships PDBs for every structure except 1tqn, and those are the only
structures whose reference output covers heteroatoms and multiple chains.

Order matters and is not incidental: MOLE draws its general-position jitter once
per atom in `structure.Atoms` order, and that order is the file's.

Out: "x y z element is_water chain seq resname bfactor atomname" per line - the
same ten columns the engine reads. Coordinates are passed through as the file's
own substrings, because the engine applies MOLE's decimal parser to them.

The atom NAME is passed rather than any derived flag: the engine holds MOLE's
backboneNames and aminoNames, so this script cannot disagree with MOLE about
what a backbone atom or a heteroatom is.
"""
import sys


def main():
    rows = 0
    with open(sys.argv[2], 'w') as o:
        for L in open(sys.argv[1]):
            if not L.startswith(('ATOM', 'HETATM')):
                continue
            name = L[12:16].strip()
            resn = L[17:20].strip()
            # Element column is optional in older PDBs; fall back to the atom
            # name's leading letters, which is what a reader has to do anyway.
            elem = L[76:78].strip() or name.lstrip('0123456789')[:1]
            chain = L[21].strip() or '?'
            # Columns 11 and 12 are the alternate-location indicator and the
            # residue insertion code, '-' when blank. MOLE's PIVOT order needs
            # both: within a residue it groups atoms by altLoc with blank first,
            # and it sorts residues by number then insertion code. Neither
            # affects the file order this script emits, which is what the
            # general-position jitter is drawn in.
            o.write('%s %s %s %s %d %s %s %s %s %s %s %s\n' % (
                L[30:38].strip(), L[38:46].strip(), L[46:54].strip(), elem,
                1 if resn == 'HOH' else 0, chain, L[22:26].strip(), resn,
                L[60:66].strip() or '0', name or '?',
                L[16].strip() or '-', L[26].strip() or '-'))
            rows += 1
    print('%d atoms -> %s' % (rows, sys.argv[2]))


main()
