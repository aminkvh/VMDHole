/* Reference dump of the cavity pipeline, for the pure-Tcl port to be checked
 * against - the stage after mole_dt_dump.
 *
 * Everything mole_build computes, at full precision (%.17g round-trips), so the
 * Tcl can be compared quantity by quantity rather than only on the cavity
 * counts at the end. The counts are a weak check: they survived several of the
 * arithmetic-form errors this port has already hit (normalize-by-multiply, the
 * snapshot in RemoveShallowVertices), because a handful of tetrahedra moving in
 * or out rarely changes how many components there are.
 *
 * Usage: mole_complex_dump ATOMS.txt OUT.txt [max_pivots]
 */
#include <stdio.h>
#include <stdlib.h>
#include "mole_complex.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M; mole_params P;
    mole_cavity *cav = NULL;
    FILE *f;
    int np, nc, nch, nvd, i, k;

    if (argc < 3) {
        fprintf(stderr, "usage: mole_complex_dump ATOMS.txt OUT.txt [max_pivots]\n");
        return 2;
    }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (argc > 3) {
        int lim = atoi(argv[3]);
        if (lim > 3 && lim < np) np = lim;
    }
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }

    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    P.min_tunnel_length = 0.0; P.weight = MOLE_W_VORONOI_SCALE;
    if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);

    f = fopen(argv[2], "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", argv[2]); return 2; }
    fprintf(f, "H %d %d %d %d %d %d\n", np, M.nt, M.n_surface, nc, nch, nvd);
    /* vdW radii: the table is part of the port and is otherwise checked only
       indirectly, through whatever elements the fixture happens to contain. */
    for (i = 0; i < np; i++) fprintf(f, "R %.17g\n", rad[i]);
    for (i = 0; i < M.nt; i++) {
        fprintf(f, "T %d %d %d %d %.17g %d", i, M.alive[i], M.boundary[i],
                M.depth[i], M.depthlen[i], M.comp[i]);
        for (k = 0; k < 3; k++) fprintf(f, " %.17g", M.center[3*i+k]);
        for (k = 0; k < 3; k++) fprintf(f, " %.17g", M.vcenter[3*i+k]);
        fprintf(f, " %.17g %.17g", M.volume[i], M.maxclear[i]);
        for (k = 0; k < 4; k++) fprintf(f, " %d", M.tv[4*i+k]);
        fprintf(f, "\n");
        fprintf(f, "E %d", i);
        for (k = 0; k < 4; k++)
            fprintf(f, " %d %.17g %.17g %.17g %.17g", M.tn[4*i+k],
                    M.eclear[4*i+k], M.elen[4*i+k], M.eweight[4*i+k], M.evweight[4*i+k]);
        fprintf(f, "\n");
    }
    for (i = 0; i < nc; i++)
        fprintf(f, "C %d %d %d %d %.17g %.17g\n", i, cav[i].count, cav[i].depth,
                cav[i].has_boundary, cav[i].volume, cav[i].depth_length);
    fclose(f);

    printf("pivots %d, tetrahedra %d, surface %d, components %d (%d channels, %d voids)\n",
           np, M.nt, M.n_surface, nc, nch, nvd);
    free(cav); mole_free(&M); dt_free(&m);
    return 0;
}
