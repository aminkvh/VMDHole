/* Reference dump of the tunnel stage, for the pure-Tcl port to be checked
 * against - the stage after mole_complex_dump.
 *
 * Dumps every INTERMEDIATE, not just the tunnels that survive: the openings,
 * the computed origins, each raw Dijkstra path, each control path, and each
 * profile's 100 sample radii and length at full precision. The tunnel list
 * alone is a weak checkpoint here for the same reason the cavity counts were
 * one stage earlier - a wrong spline or a wrong 5-nearest set still yields
 * plausible-looking tunnels of roughly the right length.
 *
 * The driver mirrors mole_main.c, which is the shipped engine, so what is
 * compared is the real call sequence rather than a test-only arrangement.
 *
 * Usage: mole_tunnel_dump ATOMS.txt OUT.txt [max_pivots] [max_similarity]
 *
 * max_similarity defaults to MOLE's 0.9. It is exposed because at 0.9 the
 * similarity filter never REMOVES anything on 1tqn or on 1BL8 - every verdict
 * is "keep" - so the sort-then-sweep removal order and the asymmetric 6-vs-2
 * sampling would go untested. 0.5 is the other end of MOLE's own clamp range
 * and does force removals.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mole_tunnel.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

typedef struct { double length; int plen, cav; mole_tunnel_profile prof; } result;

/* (Cavity.Id, Length), matching TunnelCollection.TunnelComparer and
   mole_main.c's by_cav_len - the Tcl's find_tunnels orders the same way, and
   this trace is compared against it line for line. */
static int by_cav_len(const void *a, const void *b)
{
    const result *x = a, *y = b;
    double d;
    if (x->cav != y->cav) return x->cav < y->cav ? -1 : 1;
    d = x->length - y->length;
    return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    /* zeroed, so a field added to mole_params later defaults to MOLE's own
       default here instead of to whatever was on the stack. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    int *crank = NULL;
    result *res = NULL;
    FILE *f;
    int np, nc, nch, nvd, i, k, c, nres = 0, ncap = 0;
    double max_sim = 0.9, bneck = 1.25, btol = 0.0;

    if (argc < 3) {
        fprintf(stderr, "usage: mole_tunnel_dump ATOMS.txt OUT.txt "
                        "[max_pivots] [max_sim] [bottleneck] [tolerance]\n");
        return 2;
    }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    if (argc > 3) { int lim = atoi(argv[3]); if (lim > 3 && lim < np) np = lim; }
    if (argc > 4) max_sim = atof(argv[4]);
    if (argc > 5) bneck = atof(argv[5]);
    if (argc > 6) btol = atof(argv[6]);
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }

    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    P.min_tunnel_length = 0.0; P.weight = MOLE_W_VORONOI_SCALE;
    if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);

    /* ComplexComputation.cs:436-450 numbers CHANNELS by DESCENDING VOLUME and
       assigns Cavity.Id from that; TunnelComparer then orders on the Id. It is
       NOT the order components come out of the pipeline in - on 1ERI the two
       disagree and MOLE's own ids read C1 8.99, C2 7.74, C3 11.94, C4 9.92.
       Stable, so equal volumes keep index order. */
    crank = malloc((size_t)(nc > 0 ? nc : 1) * sizeof(int));
    if (!crank) return 1;
    for (i = 0; i < nc; i++) {
        int r = 0, q;
        for (q = 0; q < nc; q++) {
            if (!cav[q].has_boundary) continue;
            if (cav[q].volume > cav[i].volume
                || (cav[q].volume == cav[i].volume && q < i)) r++;
        }
        crank[i] = r;
    }

    f = fopen(argv[2], "w");
    if (!f) { fprintf(stderr, "cannot write %s\n", argv[2]); return 2; }
    fprintf(f, "H %d %d %d %d %d\n", np, M.nt, nc, nch, nvd);

    /* Radius probes, dumped with their own coordinates so both sides evaluate
       the SAME points - hundreds of them, against the ~100 per tunnel the
       profiles give.

       What they do and do not pin, measured rather than assumed: dropping to
       the ONE nearest atom moves 151 of 600 probes, so the probes do test the
       rule. Raising it to the TEN nearest - which is what both method papers
       say, against the code - changes nothing, here or in any profile radius on
       1tqn. That is not a gap in the probes: the vdW radii of C/N/O/S span only
       ~0.3 A while the distances spread far more, so the minimum of
       (distance - vdW) is essentially never attained past the fifth atom. The
       papers' "ten" and the code's "five" are indistinguishable on protein
       data. Recorded so nobody later reads the passing test as proof of five. */
    {
        int nprobe = 0;
        for (i = 0; i < M.nt && nprobe < 600; i++) {
            if (!M.alive[i]) continue;
            fprintf(f, "V %.17g %.17g %.17g %.17g\n",
                    M.vcenter[3*i], M.vcenter[3*i+1], M.vcenter[3*i+2],
                    mole_radius_at_raw(&M.vcenter[3*i], piv, rad, np));
            fprintf(f, "V %.17g %.17g %.17g %.17g\n",
                    M.center[3*i], M.center[3*i+1], M.center[3*i+2],
                    mole_radius_at_raw(&M.center[3*i], piv, rad, np));
            nprobe += 2;
        }
    }

    for (c = 0; c < nc; c++) {
        int origins[16], nor, *open = NULL, nop, o;
        double *dist; int *prev;
        if (!(cav[c].has_boundary && cav[c].depth_length > P.min_depth_length
              && cav[c].depth > P.min_depth)) continue;
        nor = mole_auto_origins(&M, c, 10.0, 5, origins);
        fprintf(f, "G %d %d", c, nor);
        for (i = 0; i < nor; i++) fprintf(f, " %d", origins[i]);
        fprintf(f, "\n");
        if (!nor) continue;
        nop = mole_openings(&M, c, 10.0, &open);
        fprintf(f, "O %d %d", c, nop);
        for (i = 0; i < nop; i++) fprintf(f, " %d", open[i]);
        fprintf(f, "\n");
        if (!nop) { free(open); continue; }
        dist = malloc((size_t)M.nt*sizeof(double));
        prev = malloc((size_t)M.nt*sizeof(int));
        if (!dist || !prev) return 1;

        for (o = 0; o < nor; o++) {
            int base_ = nres;
            mole_dijkstra(&M, c, origins[o], dist, prev, P.weight);
            /* The whole dist/prev pair would be enormous; the paths that are
               actually followed pin the same thing and are what downstream
               reads. Reached-count and total cost catch a wrong relaxation
               that happens not to change any followed path. */
            {
                int reached = 0; double tot = 0.0;
                for (i = 0; i < M.nt; i++)
                    if (dist[i] < 1e299) { reached++; tot += dist[i]; }
                fprintf(f, "D %d %d %d %.17g\n", c, origins[o], reached, tot);
            }
            for (i = 0; i < nop; i++) {
                int v = open[i], len = 0, *p, j, ncp, *cp;
                mole_tunnel_profile pr;
                if (dist[v] >= 1e299) continue;
                for (k = v; k >= 0; k = prev[k]) len++;
                p = malloc((size_t)len*sizeof(int));
                if (!p) continue;
                for (k = v, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
                for (j = 0; j < len; j++) if (M.boundary[p[j]]) break;
                fprintf(f, "P %d %d %d %d %d", c, origins[o], v, len, j);
                for (k = 0; k < len; k++) fprintf(f, " %d", p[k]);
                fprintf(f, "\n");
                if (j < 2) { free(p); continue; }
                cp = malloc((size_t)j*sizeof(int));
                ncp = cp ? mole_control_path(&M, p, j, piv, rad, np,
                                             P.interior_threshold, cp) : 0;
                fprintf(f, "K %d %d %d %d", c, origins[o], v, ncp);
                for (k = 0; k < ncp; k++) fprintf(f, " %d", cp[k]);
                fprintf(f, "\n");
                if (ncp < 5) { free(cp); free(p); continue; }
                if (mole_profile(&M, cp, ncp, piv, rad, np, &pr) != 0) {
                    free(cp); free(p); continue;
                }
                /* The profile's own 100 samples: this is what pins the spline,
                   the 5-nearest radius query and the length integral together. */
                fprintf(f, "R %d %d %d %.17g %d", c, origins[o], v, pr.length,
                        mole_filter_bottleneck(&pr, bneck, btol, 8.0));
                for (k = 0; k < pr.sr.n; k++) fprintf(f, " %.17g", pr.sr.y[k]);
                fprintf(f, "\n");
                if (pr.length >= P.min_tunnel_length
                    && mole_filter_bottleneck(&pr, bneck, btol, 8.0)) {
                    if (nres == ncap) { ncap = ncap ? ncap*2 : 32;
                        res = realloc(res, (size_t)ncap*sizeof(*res)); }
                    res[nres].length = pr.length; res[nres].plen = j; res[nres].cav = crank[c];
                    res[nres].prof = pr; nres++;
                } else mole_profile_free(&pr);
                free(cp); free(p);
            }
            {   /* Per origin, as GetTunnels does - never across origins. */
                int cnt = nres - base_, q, w;
                if (cnt > 1) {
                    char *dd = calloc((size_t)cnt, 1);
                    mole_tunnel_profile *pf = malloc((size_t)cnt*sizeof(*pf));
                    int *pl = malloc((size_t)cnt*sizeof(int));
                    for (q = 0; q < cnt; q++) { pf[q] = res[base_+q].prof; pl[q] = res[base_+q].plen; }
                    mole_filter_similar(pf, pl, dd, cnt, max_sim);
                    fprintf(f, "S %d %d %d", c, origins[o], cnt);
                    for (q = 0; q < cnt; q++) fprintf(f, " %d", dd[q]);
                    fprintf(f, "\n");
                    w = base_;
                    for (q = 0; q < cnt; q++) {
                        if (dd[q]) mole_profile_free(&res[base_+q].prof);
                        else res[w++] = res[base_+q];
                    }
                    nres = w;
                    free(dd); free(pf); free(pl);
                }
            }
        }
        free(dist); free(prev); free(open);
    }

    qsort(res, (size_t)nres, sizeof(*res), by_cav_len);
    fprintf(f, "T %d", nres);
    for (i = 0; i < nres; i++) fprintf(f, " %.17g", res[i].length);
    fprintf(f, "\n");
    fclose(f);

    printf("pivots %d, tetrahedra %d, channels %d, tunnels %d\n", np, M.nt, nch, nres);
    for (i = 0; i < nres; i++) mole_profile_free(&res[i].prof);
    free(res); free(crank); free(cav); mole_free(&M); dt_free(&m);
    return 0;
}
