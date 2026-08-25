/* MOLE's fully automatic mode: computed origins over every cavity.
 *
 * Exists to validate the port on structures other than 1tqn without hunting for
 * a residue that happens to sit in a cavity. Reports every tunnel it finds, so
 * the result can be diffed against MOLE's own tunnels.csv.
 *
 * Usage: mole_auto_test ATOMS.txt
 *   ATOMS.txt is "x y z element is_water chain seq resname" in FILE ORDER.
 *   Set MOLE_ATOMS_EXACT=1 when the table came from MOLE's own G17 dump.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_tunnel.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];

typedef struct { double length; int cav, plen; mole_tunnel_profile prof; } result;

static int by_len(const void *a, const void *b)
{
    double d = ((const result*)a)->length - ((const result*)b)->length;
    return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    /* zeroed, so a field added to mole_params later defaults to MOLE's own
       default here instead of to whatever was on the stack. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    result *res = NULL;
    int np, nc, nch, nvd, i, k, c, nres = 0, ncap = 0, ncreated = 0, nbotrej = 0;

    if (argc < 2) { fprintf(stderr, "usage: mole_auto_test ATOMS.txt\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, NULL);
    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    P.min_tunnel_length = 0.0; P.weight = MOLE_W_VORONOI_SCALE;
    { const char *wf = getenv("MOLE_WEIGHT");
      if (wf) P.weight = (mole_weight_fn)atoi(wf);
      if (getenv("MOLE_MIN_LEN")) P.min_tunnel_length = atof(getenv("MOLE_MIN_LEN")); }
    {
        const char *ext = getenv("MOLE_TETRA");
        if (ext) {
            FILE *tf = fopen(ext, "r");
            int *tv = NULL, *tn = NULL, cap = 0, ntt = 0, a[8], q;
            if (!tf) { fprintf(stderr, "cannot open %s\n", ext); return 2; }
            while (fscanf(tf, "%d %d %d %d %d %d %d %d",
                          &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7]) == 8) {
                if (ntt == cap) { cap = cap ? cap*2 : 4096;
                    tv = realloc(tv, (size_t)cap*4*sizeof(int));
                    tn = realloc(tn, (size_t)cap*4*sizeof(int)); }
                for (q = 0; q < 4; q++) { tv[4*ntt+q] = a[q]; tn[4*ntt+q] = a[4+q]; }
                ntt++;
            }
            fclose(tf);
            fprintf(stderr, "using MOLE's triangulation: %d tetrahedra\n", ntt);
            if (mole_build_from_tetra(&M, tv, tn, ntt, piv, rad, &P) != 0) return 1;
            free(tv); free(tn);
        } else {
            if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }
            if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
        }
    }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);
    printf("atoms %d, pivots %d, tetrahedra %d, channels %d, voids %d\n",
           A.n, np, M.nt, nch, nvd);

    {   /* cavity table, ordered by descending volume as MOLE numbers them */
        int *ord = malloc((size_t)nc*sizeof(int)), nn = 0, a, b;
        for (i = 0; i < nc; i++)
            if (cav[i].has_boundary && cav[i].depth_length > P.min_depth_length
                && cav[i].depth > P.min_depth) ord[nn++] = i;
        for (a = 1; a < nn; a++) { int v = ord[a]; b = a-1;
            while (b >= 0 && cav[ord[b]].volume < cav[v].volume) { ord[b+1]=ord[b]; b--; }
            ord[b+1] = v; }
        for (a = 0; a < nn; a++) {
            int origins[16], nor2 = mole_auto_origins(&M, ord[a], 10.0, 5, origins);
            { int *op2 = NULL; int no2 = mole_openings(&M, ord[a], 10.0, &op2); free(op2);
              printf("CAV %d n=%d depth=%d origins=%d openings=%d",
                   a+1, cav[ord[a]].count, cav[ord[a]].depth, nor2, no2); }
            for (i = 0; i < nor2; i++) printf(" %d", origins[i]);
            printf("\n");
        }
        free(ord);
    }

    for (c = 0; c < nc; c++) {
        int origins[16], nor, *open = NULL, nop, o;
        double *dist; int *prev;
        if (!(cav[c].has_boundary && cav[c].depth_length > P.min_depth_length
              && cav[c].depth > P.min_depth)) continue;
        nor = mole_auto_origins(&M, c, 10.0, 5, origins);
        if (!nor) continue;
        nop = mole_openings(&M, c, 10.0, &open);
        if (!nop) { free(open); continue; }
        dist = malloc((size_t)M.nt*sizeof(double));
        prev = malloc((size_t)M.nt*sizeof(int));
        if (!dist || !prev) return 1;
        for (o = 0; o < nor; o++) {
            int base_ = nres;      /* this origin's tunnels start here */
            mole_dijkstra(&M, c, origins[o], dist, prev, P.weight);
            for (i = 0; i < nop; i++) {
                int v = open[i], len = 0, *p, j, ncp;
                int *cp;
                mole_tunnel_profile pr;
                if (dist[v] >= 1e299) continue;
                for (k = v; k >= 0; k = prev[k]) len++;
                p = malloc((size_t)len*sizeof(int));
                if (!p) continue;
                for (k = v, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
                for (j = 0; j < len; j++) if (M.boundary[p[j]]) break;
                if (j < 2) { free(p); continue; }
                ncreated++;   /* count where MOLE dumps: after TakeWhile, before controlPath */
                cp = malloc((size_t)j*sizeof(int));
                ncp = cp ? mole_control_path(&M, p, j, piv, rad, np,
                                             P.interior_threshold, cp) : 0;
                if (ncp < 5) { free(cp); free(p); continue; }
                if (mole_profile(&M, cp, ncp, piv, rad, np, &pr) != 0) {
                    free(cp); free(p); continue;
                }
                if (pr.length >= P.min_tunnel_length && mole_filter_bottleneck(&pr, 1.25, 0.0, 8.0)) {
                    if (nres == ncap) { ncap = ncap ? ncap*2 : 32;
                        res = realloc(res, (size_t)ncap*sizeof(*res)); }
                    res[nres].length = pr.length; res[nres].cav = c;
                    res[nres].plen = j; res[nres].prof = pr; nres++;
                } else { nbotrej++; mole_profile_free(&pr); }
                free(cp); free(p);
            }
            /* MOLE runs the similarity filter inside GetTunnels(origin), so it
               only ever compares tunnels that share a start point. Applying it
               across all origins compares routes from different cavities, which
               can never be similar, and lets a pivot from one origin delete
               another origin's tunnel. */
            if (getenv("MOLE_CALL_DUMP")) {
                int q;
                fprintf(stderr, "OURCALL");
                for (q = base_; q < nres; q++) fprintf(stderr, " %.2f/p%d", res[q].length, res[q].plen);
                fprintf(stderr, "\n");
            }
            {
                int cnt = nres - base_, q, w;
                if (cnt > 1) {
                    char *dd = calloc((size_t)cnt, 1);
                    mole_tunnel_profile *pf = malloc((size_t)cnt*sizeof(*pf));
                    int *pl = malloc((size_t)cnt*sizeof(int));
                    for (q = 0; q < cnt; q++) { pf[q] = res[base_+q].prof; pl[q] = res[base_+q].plen; }
                    mole_filter_similar(pf, pl, dd, cnt, 0.9);
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
    {
        result *out2 = malloc((size_t)(nres?nres:1)*sizeof(*out2));
        int kept = 0;
        for (i = 0; i < nres; i++) out2[kept++] = res[i];
        qsort(out2, (size_t)kept, sizeof(*out2), by_len);
        printf("paths created %d, rejected by bottleneck %d\n", ncreated, nbotrej);
        printf("ours generation order:");
        for (i = 0; i < nres; i++) printf(" %.2f/p%d", res[i].length, res[i].plen);
        printf("\n");
        printf("raw (pre-similarity, %d):", nres);
        { result *tmp = malloc((size_t)(nres?nres:1)*sizeof(*tmp));
          memcpy(tmp, res, (size_t)nres*sizeof(*tmp));
          qsort(tmp, (size_t)nres, sizeof(*tmp), by_len);
          for (i = 0; i < nres; i++) printf(" %.2f/p%d", tmp[i].length, tmp[i].plen);
          free(tmp); }
        printf("\n");
        /* Optional per-tunnel CSV in MOLE's own shape and density (8 points/A)
           so the profiles can be compared row by row, not just by length. */
        { const char *dir = getenv("MOLE_CSV_DIR");
          if (dir) for (i = 0; i < kept; i++) {
            char fn[512]; FILE *f2; int q, ns; double dt2, acc = 0.0, pv[3], pr2[3];
            snprintf(fn, sizeof fn, "%s/ours_tunnel_%d.csv", dir, i+1);
            f2 = fopen(fn, "w"); if (!f2) continue;
            fprintf(f2, "\"T\",\"Distance\",\"Radius\",\"X\",\"Y\",\"Z\"\n");
            ns = (int)(out2[i].prof.length * 8.0); if (ns < 1) ns = 1;
            dt2 = 1.0/ns;
            for (q = 0; q <= ns; q++) {
                double tt = dt2*q;
                pv[0]=mole_spline_eval(&out2[i].prof.sx,tt);
                pv[1]=mole_spline_eval(&out2[i].prof.sy,tt);
                pv[2]=mole_spline_eval(&out2[i].prof.sz,tt);
                if (q) acc += sqrt((pv[0]-pr2[0])*(pv[0]-pr2[0])+(pv[1]-pr2[1])*(pv[1]-pr2[1])
                                 + (pv[2]-pr2[2])*(pv[2]-pr2[2]));
                fprintf(f2, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", tt, acc,
                        mole_spline_eval(&out2[i].prof.sr,tt), pv[0], pv[1], pv[2]);
                memcpy(pr2, pv, sizeof pr2);
            }
            fclose(f2);
          } }
        printf("tunnels %d (from %d before the similarity filter)\n", kept, nres);
        printf("lengths:");
        for (i = 0; i < kept; i++) printf(" %.2f", out2[i].length);
        printf("\n");
        free(out2);
    }
    for (i = 0; i < nres; i++) mole_profile_free(&res[i].prof);
    free(res); free(cav); mole_free(&M);
    return 0;
}
