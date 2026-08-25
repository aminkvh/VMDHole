/* Does our Delaunay reproduce MOLE's triangulation?
 *
 * MOLE reports 25189 tetrahedra / 50283 edges on 1tqn. Ours comes out at 25181
 * and is the EXACT Delaunay triangulation (dt_verify: 0 violations); MOLE's is
 * not, because DHTriangulation tests a raw double determinant sign.
 *
 * Usage: mole_tri_test ATOMS.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include "mole_complex.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dt_mesh m;
    int np, nt, ne = 0, i, k;

    if (argc < 2) { fprintf(stderr, "usage: mole_tri_test ATOMS.txt\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }

    nt = dt_count_finite(&m);
    for (i = 0; i < m.nt; i++) {
        if (m.t[i].dead || !dt_is_finite(&m, i)) continue;
        for (k = 0; k < 4; k++) {
            int nb = m.t[i].nb[k];
            if (nb >= 0 && !m.t[nb].dead && dt_is_finite(&m, nb) && nb > i) ne++;
        }
    }
    printf("atoms read      %d\n", A.n);
    printf("pivots          %d          (expect 3809)\n", np);
    printf("tetrahedra      %d          (MOLE 25189, theirs is not Delaunay)\n", nt);
    printf("edges           %d          (MOLE 50283)\n", ne);
    dt_free(&m);
    return 0;
}
