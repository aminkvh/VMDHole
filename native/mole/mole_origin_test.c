/* Checkpoint: does the ported origin snap land where MOLE's does?
 *
 * Reference (test.xml, origin = residues A308 + A309):
 *   Id User1212 A308 A
 *   Cavity 1, Delta = 1.17ang, Near {ARG 212 A, GLU 308 A}
 *
 * A sharp check - it pins a cavity, a distance, and a residue pair at once. Two
 * details in Computation.cs decide whether it can pass at all:
 *   - the vertex is SELECTED by distance to VoronoiCenter, but the Delta that
 *     gets printed is measured to Center. Different points.
 *   - the requested point is the centroid of every atom of the listed residues,
 *     taken AFTER the jitter.
 *
 * Usage: mole_origin_test ATOMS.txt CHAIN:SEQ [CHAIN:SEQ ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_rng.h"
#include "mole_complex.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];
static int    pivsrc[MOLE_MAXATOM];

typedef struct { int idx; double vol; } cav_order;

static int by_vol_desc(const void *a, const void *b)
{
    double d = ((const cav_order*)b)->vol - ((const cav_order*)a)->vol;
    return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    /* zeroed, so a field added to mole_params later defaults to MOLE's own
       default here instead of to whatever was on the stack. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    cav_order *ord = NULL;
    double c[3] = {0,0,0};
    int np, i, k, nc, nch, nvd, nsel = 0, nchan = 0;

    if (argc < 3) { fprintf(stderr, "usage: mole_origin_test ATOMS.txt CHAIN:SEQ ...\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, pivsrc);

    /* Origin point: centroid of every atom of the requested residues. */
    for (k = 2; k < argc; k++) {
        char cc[8]; int ss;
        if (sscanf(argv[k], "%7[^:]:%d", cc, &ss) != 2) continue;
        for (i = 0; i < A.n; i++)
            if (A.seq[i] == ss && !strcmp(A.chain[i], cc)) {
                c[0] += A.xyz[3*i]; c[1] += A.xyz[3*i+1]; c[2] += A.xyz[3*i+2]; nsel++;
            }
    }
    if (!nsel) { fprintf(stderr, "no atoms matched\n"); return 1; }
    for (k = 0; k < 3; k++) c[k] /= nsel;

    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }
    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);

    /* Cavity ids follow descending volume among the channels only. */
    ord = malloc((size_t)nc * sizeof(*ord));
    for (i = 0; i < nc; i++)
        if (cav[i].has_boundary && cav[i].depth_length > P.min_depth_length
            && cav[i].depth > P.min_depth) {
            ord[nchan].idx = i; ord[nchan].vol = cav[i].volume; nchan++;
        }
    qsort(ord, (size_t)nchan, sizeof(*ord), by_vol_desc);

    printf("origin point   %.3f %.3f %.3f   (centroid of %d atoms)\n", c[0], c[1], c[2], nsel);
    for (k = 0; k < nchan; k++) {
        int best = -1, t;
        double bd = 1e300;
        for (t = 0; t < M.nt; t++) {
            double dx, dy, dz, d2;
            if (!M.alive[t] || M.comp[t] != ord[k].idx) continue;
            if (M.depth[t] < 5) continue;              /* GetOrigin's own filter */
            dx = M.vcenter[3*t]-c[0]; dy = M.vcenter[3*t+1]-c[1]; dz = M.vcenter[3*t+2]-c[2];
            d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < bd) { bd = d2; best = t; }
        }
        if (best < 0 || bd > 5.0*5.0) continue;        /* OriginRadius */
        {
            double dx = M.center[3*best]-c[0], dy = M.center[3*best+1]-c[1],
                   dz = M.center[3*best+2]-c[2];
            char seen[4][24]; int ns = 0, j;
            printf("  => Cavity %d, Delta = %.2fang, Near {", k+1,
                   sqrt(dx*dx+dy*dy+dz*dz));
            for (i = 0; i < 4; i++) {
                int a = pivsrc[M.tv[4*best+i]];
                char lab[24];
                snprintf(lab, sizeof lab, "%s %d %s", A.resn[a], A.seq[a], A.chain[a]);
                for (j = 0; j < ns; j++) if (!strcmp(seen[j], lab)) break;
                if (j == ns) { snprintf(seen[ns], sizeof seen[ns], "%s", lab); ns++; }
            }
            for (j = 0; j < ns; j++) printf("%s%s", j ? ", " : "", seen[j]);
            printf("}\n");
        }
    }
    printf("\nMOLE:  Cavity 1, Delta = 1.17ang, Near {ARG 212 A, GLU 308 A}\n");
    free(ord); free(cav); mole_free(&M); dt_free(&m);
    return 0;
}
