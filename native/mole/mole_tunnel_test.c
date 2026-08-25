/* Checkpoint: does the ported tunnel stage find MOLE's 3 tunnels on 1tqn?
 *
 * Reference (test.xml, origin = residues A308 + A309): 3 tunnels, lengths
 * 13.82 / 29.89 / 34.98.
 *
 * Usage: mole_tunnel_test ATOMS.txt CHAIN:SEQ [CHAIN:SEQ ...]
 *
 * With MOLE_TETRA=<file> the triangulation is read from MOLE's own dump
 * ("v0 v1 v2 v3 n0 n1 n2 n3" per tetrahedron) instead of being computed here.
 * That isolates the port from the 0.03% triangulation difference: if the output
 * then matches MOLE's, everything downstream of the triangulation is exact.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_rng.h"
#include "mole_tunnel.h"

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];
static int    pivsrc[MOLE_MAXATOM];
#define MAXTUN  512

typedef struct { int *path, np, dead; mole_tunnel_profile prof; } tun;

/* Sample a profile's centreline at MOLE's own density for the similarity test:
   GetProfile(6), i.e. 6 points per Angstrom of tunnel length. */
static int sample_line(const mole_tunnel_profile *p, double density, double **out)
{
    int n = (int)(p->length * density), i;
    double dt;
    if (n < 1) n = 1;
    dt = 1.0 / n;
    *out = malloc((size_t)(n + 1) * 3 * sizeof(double));
    if (!*out) return 0;
    for (i = 0; i <= n; i++) {
        double t = dt * i;
        (*out)[3*i+0] = mole_spline_eval(&p->sx, t);
        (*out)[3*i+1] = mole_spline_eval(&p->sy, t);
        (*out)[3*i+2] = mole_spline_eval(&p->sz, t);
    }
    return n + 1;
}

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    /* zeroed, so a field added to mole_params later defaults to MOLE's own
       default here instead of to whatever was on the stack. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    static tun T[MAXTUN];
    double c[3] = {0,0,0}, *dist, bestvol = -1;
    int *prev, *open = NULL;
    int np, i, j, k, nc, nch, nvd, nsel = 0;
    int cav1 = -1, origin = -1, nop, nt = 0, kept = 0;

    if (argc < 3) { fprintf(stderr, "usage: mole_tunnel_test ATOMS.txt CHAIN:SEQ ...\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, pivsrc);
    for (k = 2; k < argc; k++) {
        char cc[8]; int ss;
        if (sscanf(argv[k], "%7[^:]:%d", cc, &ss) != 2) continue;
        for (i = 0; i < A.n; i++)
            if (A.seq[i] == ss && !strcmp(A.chain[i], cc)) {
                c[0]+=A.xyz[3*i]; c[1]+=A.xyz[3*i+1]; c[2]+=A.xyz[3*i+2]; nsel++;
            }
    }
    if (!nsel) { fprintf(stderr, "no atoms matched\n"); return 1; }
    for (k = 0; k < 3; k++) c[k] /= nsel;

    P.probe_radius=3.0; P.interior_threshold=1.25; P.min_depth=8; P.min_depth_length=5.0;
    P.min_tunnel_length = 0.0; P.weight = MOLE_W_VORONOI_SCALE;
    { const char *wf = getenv("MOLE_WEIGHT");
      if (wf) P.weight = (mole_weight_fn)atoi(wf);
      if (getenv("MOLE_MIN_LEN")) P.min_tunnel_length = atof(getenv("MOLE_MIN_LEN")); }
    {
        const char *ext = getenv("MOLE_TETRA");
        if (ext) {
            FILE *tf = fopen(ext, "r");
            int *tv = NULL, *tn = NULL, cap = 0, ntt = 0, a[8];
            if (!tf) { fprintf(stderr, "cannot open %s\n", ext); return 2; }
            while (fscanf(tf, "%d %d %d %d %d %d %d %d",
                          &a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7]) == 8) {
                if (ntt == cap) { cap = cap ? cap*2 : 4096;
                    tv = realloc(tv, (size_t)cap*4*sizeof(int));
                    tn = realloc(tn, (size_t)cap*4*sizeof(int)); }
                for (k = 0; k < 4; k++) { tv[4*ntt+k] = a[k]; tn[4*ntt+k] = a[4+k]; }
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

    for (i = 0; i < nc; i++)
        if (cav[i].has_boundary && cav[i].depth_length > P.min_depth_length
            && cav[i].depth > P.min_depth && cav[i].volume > bestvol) {
            bestvol = cav[i].volume; cav1 = i;
        }
    {
        double bd = 1e300;
        for (i = 0; i < M.nt; i++) {
            double dx, dy, dz, d2;
            if (!M.alive[i] || M.comp[i] != cav1 || M.depth[i] < 5) continue;
            dx=M.vcenter[3*i]-c[0]; dy=M.vcenter[3*i+1]-c[1]; dz=M.vcenter[3*i+2]-c[2];
            d2 = dx*dx+dy*dy+dz*dz;
            if (d2 < bd) { bd = d2; origin = i; }
        }
    }
    if (origin < 0) { fprintf(stderr, "no origin\n"); return 1; }

    nop = mole_openings(&M, cav1, 10.0, &open);
    /* MOLE picks its opening pivots with Rx's MaxBy over a HashSet<Tetrahedron>
       that has no GetHashCode override, so ties are broken by .NET object-identity
       hashing - allocation order, not geometry. MOLE_OPENINGS lets their choice
       be substituted to test whether that tie is the only thing left. */
    {
        const char *ov = getenv("MOLE_OPENINGS");
        if (ov) {
            int *o2 = malloc(64*sizeof(int)); int n2 = 0;
            char buf[512]; strncpy(buf, ov, sizeof buf - 1); buf[sizeof buf - 1] = 0;
            { char *tok = strtok(buf, ", ");
              while (tok && n2 < 64) { o2[n2++] = atoi(tok); tok = strtok(NULL, ", "); } }
            free(open); open = o2; nop = n2;
            fprintf(stderr, "forced %d openings\n", nop);
        }
    }
    dist = malloc((size_t)M.nt*sizeof(double));
    prev = malloc((size_t)M.nt*sizeof(int));
    if (!dist || !prev) return 1;
    mole_dijkstra(&M, cav1, origin, dist, prev, P.weight);

    for (i = 0; i < nop && nt < MAXTUN; i++) {
        int v = open[i], len = 0, *p;
        if (dist[v] >= 1e299) continue;
        for (k = v; k >= 0; k = prev[k]) len++;
        p = malloc((size_t)len*sizeof(int));
        if (!p) continue;
        for (k = v, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
        /* Tunnel.Create: path = path.TakeWhile(p => !p.IsBoundary) */
        for (j = 0; j < len; j++) if (M.boundary[p[j]]) break;
        if (j < 2) { free(p); continue; }
        /* The spline is fitted to the CONTROL path, not the Dijkstra path. */
        {
            int *cp = malloc((size_t)j * sizeof(int));
            int ncp = cp ? mole_control_path(&M, p, j, piv, rad, np,
                                             P.interior_threshold, cp) : 0;
            if (ncp < 5) { free(cp); free(p); continue; }
            T[nt].path = p; T[nt].np = j; T[nt].dead = 0;
            if (mole_profile(&M, cp, ncp, piv, rad, np, &T[nt].prof) != 0) {
                free(cp); free(p); continue;
            }
            free(cp);
        }
        if (T[nt].prof.length < P.min_tunnel_length
            || !mole_filter_bottleneck(&T[nt].prof, 1.25, 0.0, 8.0)) { free(p); continue; }
        nt++;
    }

    /* Complex.RemoveLonger: for a pair, take the one with the SHORTER path and
       ask what fraction of its 6/A points lie within 1 A of the other. Above
       MaxTunnelSimilarity the longer one goes. */
    for (i = 0; i < nt; i++)
        for (j = i+1; j < nt; j++) {
            int a = i, b = j, na, nb, q, s, hit = 0;
            double *la, *lb;
            if (T[j].np < T[i].np) { a = j; b = i; }
            if (T[a].dead || T[b].dead) continue;
            na = sample_line(&T[a].prof, 6.0, &la);
            nb = sample_line(&T[b].prof, 6.0, &lb);
            if (!na || !nb) { free(la); free(lb); continue; }
            for (q = 0; q < na; q++) {
                double best = 1e300;
                for (s = 0; s < nb; s++) {
                    double dx=la[3*q]-lb[3*s], dy=la[3*q+1]-lb[3*s+1], dz=la[3*q+2]-lb[3*s+2];
                    double d2 = dx*dx+dy*dy+dz*dz;
                    if (d2 < best) best = d2;
                }
                if (best <= 1.0) hit++;
            }
            if ((double)hit/na > 0.9) T[b].dead = 1;
            free(la); free(lb);
        }
    for (i = 0; i < nt; i++) if (!T[i].dead) kept++;

    /* Graph dump, for asking whether MOLE's route even EXISTS in our graph.
       If it does and is cheaper than ours, our search is wrong; if it does and
       is dearer, our search is right and the weights or the graph differ; if it
       does not, the difference is upstream of the search entirely. */
    if (getenv("MOLE_DUMP_GRAPH")) {
        FILE *g = fopen("/tmp/ours_graph.txt", "w");
        if (g) {
            for (i = 0; i < M.nt; i++) {
                if (!M.alive[i] || M.comp[i] != cav1) continue;
                fprintf(g, "V %d %.6f %.6f %.6f %d %d\n", i,
                        M.vcenter[3*i], M.vcenter[3*i+1], M.vcenter[3*i+2],
                        M.boundary[i], M.depth[i]);
            }
            for (i = 0; i < M.nt; i++) {
                if (!M.alive[i] || M.comp[i] != cav1) continue;
                for (k = 0; k < 4; k++) {
                    int v = M.tn[4*i+k];
                    if (v < 0 || !M.alive[v] || M.comp[v] != cav1 || v < i) continue;
                    fprintf(g, "E %d %d %.17g %.17g\n", i, v, M.eweight[4*i+k], M.eclear[4*i+k]);
                }
            }
            fprintf(g, "K %d", nop);
            for (i = 0; i < nop; i++) fprintf(g, " %d", open[i]);
            fprintf(g, "\n");
            for (i = 0; i < nt; i++) {
                if (T[i].dead) continue;
                fprintf(g, "P");
                for (j = 0; j < T[i].np; j++) fprintf(g, " %d", T[i].path[j]);
                fprintf(g, "\n");
            }
            fprintf(g, "O %d\n", origin);
            fclose(g);
        }
    }

    printf("cavity 1 openings   %d\n", nop);
    printf("raw paths           %d   (after boundary truncation + bottleneck filter)\n", nt);
    printf("tunnels             %d   (MOLE 3)   %s\n", kept, kept == 3 ? "MATCH" : "MISMATCH");
    /* Dump our profiles in MOLE's own CSV shape and density (8 points/A, which
       is what their 111 rows over 13.817 A works out to) so the two can be
       compared row for row rather than by summary numbers. */
    {
        int id = 0;
        for (i = 0; i < nt; i++) {
            char fn[256]; FILE *o; int q, ns;
            double dt2, acc = 0.0, pxv[3], prv[3];
            if (T[i].dead) continue;
            snprintf(fn, sizeof fn, "/tmp/ours_tunnel_%d.csv", ++id);
            o = fopen(fn, "w");
            if (!o) continue;
            fprintf(o, "\"T\",\"Distance\",\"Radius\",\"X\",\"Y\",\"Z\"\n");
            ns = (int)(T[i].prof.length * 8.0); if (ns < 1) ns = 1;
            dt2 = 1.0/ns;
            for (q = 0; q <= ns; q++) {
                double tt = dt2*q;
                pxv[0]=mole_spline_eval(&T[i].prof.sx,tt);
                pxv[1]=mole_spline_eval(&T[i].prof.sy,tt);
                pxv[2]=mole_spline_eval(&T[i].prof.sz,tt);
                if (q) acc += sqrt((pxv[0]-prv[0])*(pxv[0]-prv[0])
                                 + (pxv[1]-prv[1])*(pxv[1]-prv[1])
                                 + (pxv[2]-prv[2])*(pxv[2]-prv[2]));
                fprintf(o, "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", tt, acc,
                        mole_spline_eval(&T[i].prof.sr,tt), pxv[0], pxv[1], pxv[2]);
                memcpy(prv, pxv, sizeof prv);
            }
            fclose(o);
        }
    }

    printf("\n  ours lengths: ");
    for (i = 0; i < nt; i++) if (!T[i].dead) printf("%.2f ", T[i].prof.length);
    printf("\n  MOLE lengths: 13.82 29.89 34.98\n");
    return kept == 3 ? 0 : 1;
}
