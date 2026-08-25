/* Chain-collation regression for two sites a fixture-only check cannot
 * distinguish: ResidueFlow (mole_build_flow) and HetResidues
 * (mole_sort_residues_by_chain). Both route through mole_chain_cmp; this pins
 * them with residue tables built to actually distinguish it, rather than
 * fixtures that pass under either rule.
 *
 * Case 1 - ResidueFlow. MOLE's first layer takes the five nearest atoms with
 * no chain restriction, so oligomeric pores commonly have multi-chain first
 * layers: 1BL8's own reference has THR 107 A, THR 107 C, THR 107 D as its
 * first ResidueFlow entries (tunnels.xml). Relabel chain A -> a: three
 * DIFFERENT letters (a, c, d), so InvariantCulture's alphabetical-by-letter
 * rule gives a, C, D - where ordinal strcmp gives C, D, a.
 *
 * Case 2 - HetResidues. 1ERI tunnel 3 carries DA 6 B, DA 7 B (tunnels.xml).
 * Relabel the DA 6 chain to lowercase b: chain 'b' and chain 'B' are the SAME
 * letter, so this is the collation tiebreak (lowercase first) rather than
 * letter order - a different half of mole_chain_cmp than case 1 exercises.
 * DA6(b) must sort strictly before DA7(B) regardless of residue number.
 */
#include <stdio.h>
#include <string.h>
#include "mole_lining.h"

static void set_res(mole_residues *R, int i, const char *name, const char *chain, int seq)
{
    strcpy(R->name[i], name);
    strcpy(R->chain[i], chain);
    R->seq[i] = seq;
}

static int case_residue_flow(void)
{
    mole_residues R;
    mole_layer L[1];
    mole_flow_entry flow[3];
    int nflow, bad = 0;
    R.n = 3;
    /* File/input order D, C, a - deliberately not already sorted. */
    set_res(&R, 0, "THR", "D", 107);
    set_res(&R, 1, "THR", "C", 107);
    set_res(&R, 2, "THR", "a", 107);

    memset(&L[0], 0, sizeof L[0]);
    L[0].nres = 3;
    L[0].res[0] = 0; L[0].isbb[0] = 0;
    L[0].res[1] = 1; L[0].isbb[1] = 0;
    L[0].res[2] = 2; L[0].isbb[2] = 0;

    nflow = mole_build_flow(L, 1, &R, flow);
    printf("  ResidueFlow: got %d entries, order", nflow);
    { int i; for (i = 0; i < nflow; i++) printf(" %s%d%s", R.name[flow[i].res],
              R.seq[flow[i].res], R.chain[flow[i].res]); }
    printf("\n");
    if (nflow != 3 || strcmp(R.chain[flow[0].res], "a") || strcmp(R.chain[flow[1].res], "C")
        || strcmp(R.chain[flow[2].res], "D")) {
        printf("  FAIL: want order a, C, D (InvariantCulture) - THR 107 A/C/D relabelled a/C/D\n");
        bad = 1;
    }
    return bad;
}

static int case_het_residues(void)
{
    mole_residues R;
    int out[2], bad = 0;
    R.n = 2;
    set_res(&R, 0, "DA", "b", 6);
    set_res(&R, 1, "DA", "B", 7);
    /* Input order DA6(b), DA7(B) - already the collated order, so a correct
       sort leaves it unchanged; an ordinal sort (strcmp 'B' < 'b') swaps it. */
    out[0] = 0; out[1] = 1;

    mole_sort_residues_by_chain(out, 2, &R);
    printf("  HetResidues: got order DA%d%s, DA%d%s\n",
           R.seq[out[0]], R.chain[out[0]], R.seq[out[1]], R.chain[out[1]]);
    if (out[0] != 0 || out[1] != 1) {
        printf("  FAIL: want DA6 b, DA7 B unchanged (lowercase sorts before uppercase)\n");
        bad = 1;
    }
    return bad;
}

int main(void)
{
    int bad = 0;
    bad |= case_residue_flow();
    bad |= case_het_residues();
    if (bad) { printf("  FAIL\n"); return 1; }
    printf("  PASS\n");
    return 0;
}
