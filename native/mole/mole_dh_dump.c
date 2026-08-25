/* Dump the DH triangulation in MOLE's own MOLE_TETRA_DUMP format, so the port
   can be compared against MOLE cell for cell and SLOT for slot before anything
   is wired into the engine.

   Centroids are printed raw; the comparator normalises "-0.0000", which C's
   printf produces for a small negative where .NET's formatter gives "0.0000".

   Usage: mole_dh_dump ATOMS.txt > cells.txt */
#include <stdio.h>
#include <stdlib.h>
#include "mole_complex.h"
#include "mole_dh.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dh_mesh m;
    int np, t, k;

    if (argc < 2) { fprintf(stderr, "usage: %s ATOMS.txt\n", argv[0]); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (dh_build(&m, piv, np) != 0) { fprintf(stderr, "dh_build failed\n"); return 1; }

    fprintf(stderr, "dh: %d pivots, %d finite cells\n", np, m.nt);
    for (t = 0; t < m.nt; t++) {
        double cx = 0, cy = 0, cz = 0;
        for (k = 0; k < 4; k++) {
            int a = m.tv[4*t+k];
            cx += piv[3*a]; cy += piv[3*a+1]; cz += piv[3*a+2];
        }
        {   /* MOLE's third field: how many of the four neighbours are absent
               (its infinite cells are dropped, so its Adjacency entry is null). */
            int nnull = 0;
            for (k = 0; k < 4; k++) if (m.tn[4*t+k] < 0) nnull++;
            printf("%.4f %.4f %.4f %d", cx/4, cy/4, cz/4, nnull);
        }
        for (k = 0; k < 4; k++) {
            int a = m.tv[4*t+k];
            printf("%c%.3f,%.3f,%.3f", k ? ';' : ' ',
                   piv[3*a], piv[3*a+1], piv[3*a+2]);
        }
        putchar('\n');
    }
    dh_mesh_free(&m);
    return 0;
}
