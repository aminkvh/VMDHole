/* Checkpoint: does the ported pipeline find MOLE's cavities on 1tqn?
 * Reference (test.xml): 9 channels, 6 voids.
 * Usage: mole_cavity_test ATOMS.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include "mole_complex.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    /* zeroed, so a field added to mole_params later defaults to MOLE's own
       default here instead of to whatever was on the stack. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    int np, nc, nch, nvd;

    if (argc < 2) { fprintf(stderr, "usage: mole_cavity_test ATOMS.txt\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }

    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);

    printf("pivots            %d\n", np);
    printf("tetrahedra        %d          (MOLE 25189)\n", M.nt);
    printf("surface cavity    %d tetrahedra after the probe peel\n", M.n_surface);
    printf("components        %d\n", nc);
    printf("channels          %d          (MOLE 9)   %s\n", nch, nch == 9 ? "MATCH" : "MISMATCH");
    printf("voids             %d          (MOLE 6)   %s\n", nvd, nvd == 6 ? "MATCH" : "MISMATCH");
    free(cav); mole_free(&M); dt_free(&m);
    return (nch == 9 && nvd == 6) ? 0 : 1;
}
