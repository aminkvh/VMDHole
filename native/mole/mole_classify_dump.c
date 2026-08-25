/* Dump mole_is_backbone_name / mole_is_amino_name for every name the Tcl side
   checks, so the two copies of MOLE's sets can be diffed.
 *
 * PdbEx.backboneNames and PdbResidue.aminoNames are the only part of the port
 * with no ground-truth coverage: every reference structure is protein, single
 * model, no hydrogens and no nucleic acid, so the H entry and the whole nucleic
 * half are never exercised by a MOLE comparison. Two hand-kept copies of a set
 * drift; this is what stops them.
 *
 * Built and run by tests/test_mole_tcl_port.sh, not by build.sh - it is a test
 * fixture, not part of the engine. */
#include <stdio.h>
#include "mole_complex.h"
int main(void)
{
    static const char *N[] = {
        "C","N","O","H","CA","P","O1P","O2P","OP1","OP2","O5'","C5'","C4'",
        "O4'","C1'","C2'","C3'","O3'","O2'","CB","CG","CG1","CG2","OG","OG1",
        "CD","CD1","NE","NZ","OH","SD","FE","ZN","MG","HOH","XX","ca","o5'",
        "ALA","ARG","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
        "MET","PHE","PRO","SER","THR","TRP","TYR","VAL","ASN","HEM","TIP3",
        "WAT","SOL","DA","DC","DG","DT","U","MSE","UNK","his","hem", NULL };
    int i;
    for (i = 0; N[i]; i++) printf("BB %s %d\n", N[i], mole_is_backbone_name(N[i]));
    for (i = 0; N[i]; i++) printf("AA %s %d\n", N[i], mole_is_amino_name(N[i]));
    return 0;
}
