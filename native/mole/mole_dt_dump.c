/* Reference dump of the Delaunay build, for the pure-Tcl port to be checked
 * against.
 *
 * The check is deliberately split in two, because a single combined comparison
 * cannot say which half broke:
 *
 *   POINTS  the quantised integer coordinates, super-tetrahedron corners
 *           included. Tests jitter + quantisation only.
 *   TETRA   every tetrahedron SLOT, dead ones included, with its neighbours,
 *           plus the free list. Tests the Bowyer-Watson logic only, since the
 *           Tcl side is fed the POINTS integers rather than re-deriving them.
 *
 * Slots and indices, not the live set: tetrahedron indices are load-bearing
 * downstream (component iteration is index-ordered, and a MaxClearance tie is
 * settled by lowest index), so matching the live set as a SET would pass while
 * still producing different tunnels.
 *
 * Usage: mole_dt_dump ATOMS.txt POINTS.txt TETRA.txt [max_pivots]
 *   max_pivots truncates the point set, so parity can be established and timed
 *   on a few hundred points before a whole protein is attempted.
 */
#include <stdio.h>
#include <stdlib.h>
#include "mole_complex.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dt_mesh m;
    FILE *fp, *ft;
    int np, i, k;

    if (argc < 4) {
        fprintf(stderr, "usage: mole_dt_dump ATOMS.txt POINTS.txt TETRA.txt [max_pivots]\n");
        return 2;
    }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (argc > 4) {
        int lim = atoi(argv[4]);
        if (lim > 3 && lim < np) np = lim;
    }
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }

    fp = fopen(argv[2], "w");
    if (!fp) { fprintf(stderr, "cannot write %s\n", argv[2]); return 2; }
    fprintf(fp, "%d\n", m.npt);
    for (i = 0; i < m.npt + 4; i++)
        fprintf(fp, "%ld %ld %ld\n", m.p[3*i], m.p[3*i+1], m.p[3*i+2]);
    fclose(fp);

    ft = fopen(argv[3], "w");
    if (!ft) { fprintf(stderr, "cannot write %s\n", argv[3]); return 2; }
    fprintf(ft, "%d %d %d\n", m.nt, m.nfree, m.last);
    for (i = 0; i < m.nt; i++) {
        fprintf(ft, "%d", m.t[i].dead);
        for (k = 0; k < 4; k++) fprintf(ft, " %d", m.t[i].v[k]);
        for (k = 0; k < 4; k++) fprintf(ft, " %d", m.t[i].nb[k]);
        fprintf(ft, "\n");
    }
    for (i = 0; i < m.nfree; i++) fprintf(ft, "%d\n", m.freet[i]);
    fclose(ft);

    printf("pivots %d, tetrahedra %d slots (%d live, %d finite), free %d\n",
           np, m.nt, dt_count_live(&m), dt_count_finite(&m), m.nfree);
    dt_free(&m);
    return 0;
}
