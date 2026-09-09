/* conn_lobes.c - fast replacement for the plugin's two pure-Tcl hot loops in
 * the Connolly lateral-opening ("lobe") coloring path: classifying every
 * .sph dot as pore-or-lateral and clustering the lateral dots into lobes
 * (proc ::VMDPathFinder::_conn_classify_sph + _conn_frame_lobes), and assigning
 * each triangle of the unified surface mesh to the region its nearest dot
 * belongs to (proc ::VMDPathFinder::_split_conn_mesh_by_region). Measured on a
 * real 108,650-dot Connolly frame: classify+cluster 1.34s, split 2.12s in
 * Tcl - both are the SAME arithmetic this file does, line for line, just
 * compiled instead of interpreted with per-token string operations.
 *
 * Two subcommands, matching the two call sites:
 *
 *   conn_lobes classify SPH CX CY CZ VX VY VZ MARGIN [F1X F1Y F1Z F2X F2Y F2Z]
 *     Prints PORE/KEEP/LATERAL/ESCRANGE/LOBE/MARKED blocks (see usage()) that a Tcl
 *     wrapper turns back into the exact dict shapes _conn_classify_sph and
 *     _conn_frame_lobes already return, so every downstream Tcl caller
 *     (site-table discovery, per-frame region building, two-tone coloring)
 *     is unchanged. MARKED lists the records a LAST-REC-END followed in the
 *     source (HOLE's clip-sphere marker); writers put it back after them.
 *
 *   conn_lobes split UNION_PLOT --region NAME LABELS_SPH OUT_PLOT [...]
 *     LABELS_SPH is a plain .sph with one ATOM/HETATM line per dot in that
 *     region (only x,y,z read) - the same "lines" list the Tcl split
 *     function classifies against, not the mesh-building rsph (which adds
 *     centreline points for the pore region that must NOT enter the
 *     classification grid). Writes each OUT_PLOT with "draw delete all"
 *     plus every triangle line whose centroid's nearest labelled dot is in
 *     that region, byte-identical to the corresponding input line.
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hole_io.h"
#include <math.h>
#include <limits.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void *xmalloc(size_t n) { void *p = malloc(n); if (!p) { fprintf(stderr, "out of memory\n"); exit(1); } return p; }
static void *xrealloc(void *p, size_t n) { void *q = realloc(p, n); if (!q) { fprintf(stderr, "out of memory\n"); exit(1); } return q; }

static char *xstrdup(const char *s) {
    size_t n = strlen(s);
    char *p = xmalloc(n + 1);
    memcpy(p, s, n + 1);
    return p;
}

/* ------------------------------------------------------------- PDB fields */

/* Tcl's [string range $line A B] is inclusive; a field spanning columns A..B
   (0-indexed) is B-A+1 characters. Out-of-range or short lines read as
   blank, matching Tcl's own out-of-range behaviour (never errors). */
static void field(const char *line, size_t len, int a, int b, char *out, int outsz) {
    int n = b - a + 1, k = 0;
    for (int i = 0; i < n && k < outsz - 1; i++) {
        int p = a + i;
        out[k++] = (p >= 0 && (size_t)p < len) ? line[p] : ' ';
    }
    out[k] = 0;
    while (k > 0 && (out[k-1] == ' ' || out[k-1] == '\t' || out[k-1] == '\r')) out[--k] = 0;
    int s = 0; while (out[s] == ' ') s++;
    if (s) memmove(out, out + s, strlen(out + s) + 1);
}

static int is_double(const char *s, double *v) {
    if (!*s) return 0;
    char *end;
    *v = strtod(s, &end);
    while (*end == ' ') end++;
    return *end == 0;
}

static int is_atom_line(const char *line) {
    return !strncmp(line, "ATOM  ", 6) || !strncmp(line, "HETATM", 6);
}

/* -------------------------------------------------------------- dot list */

typedef struct { double t, x, y, z; char *line; } Dot;
typedef struct { double t, x, y, z, r; } Cen;

static Dot *g_dots = NULL; int g_ndots = 0, g_dotcap = 0;
/* LAST-REC-END: which stored record it followed (dot index >= 0, keep index
   encoded as -(k+1)), so writers can put the marker back. */
static char *g_dotmark = NULL; static char *g_keepmark = NULL;
static char **g_keep = NULL; int g_nkeep = 0, g_keepcap = 0;
static Cen *g_cen = NULL; int g_ncen = 0, g_cencap = 0;
static double *g_esct = NULL; int g_nesc = 0, g_esccap = 0;

static void push_dot(double t, double x, double y, double z, const char *line) {
    if (g_ndots >= g_dotcap) {
        int old = g_dotcap; g_dotcap = g_dotcap ? g_dotcap * 2 : 4096;
        g_dots = xrealloc(g_dots, g_dotcap * sizeof(Dot));
        g_dotmark = xrealloc(g_dotmark, g_dotcap); memset(g_dotmark + old, 0, g_dotcap - old);
    }
    g_dots[g_ndots].t = t; g_dots[g_ndots].x = x; g_dots[g_ndots].y = y; g_dots[g_ndots].z = z;
    g_dots[g_ndots].line = xstrdup(line);
    g_ndots++;
}
static void mark_last(int which, int idx) {
    if (which == 0) g_dotmark[idx] = 1; else g_keepmark[idx] = 1;
}
static void push_keep(const char *line) {
    if (g_nkeep >= g_keepcap) {
        int old = g_keepcap; g_keepcap = g_keepcap ? g_keepcap * 2 : 256;
        g_keep = xrealloc(g_keep, g_keepcap * sizeof(char *));
        g_keepmark = xrealloc(g_keepmark, g_keepcap); memset(g_keepmark + old, 0, g_keepcap - old);
    }
    g_keep[g_nkeep++] = xstrdup(line);
}
static void push_cen(double t, double x, double y, double z, double r) {
    if (g_ncen >= g_cencap) { g_cencap = g_cencap ? g_cencap * 2 : 1024; g_cen = xrealloc(g_cen, g_cencap * sizeof(Cen)); }
    g_cen[g_ncen].t = t; g_cen[g_ncen].x = x; g_cen[g_ncen].y = y; g_cen[g_ncen].z = z; g_cen[g_ncen].r = r;
    g_ncen++;
}
static void push_esc(double t) {
    if (g_nesc >= g_esccap) { g_esccap = g_esccap ? g_esccap * 2 : 256; g_esct = xrealloc(g_esct, g_esccap * sizeof(double)); }
    g_esct[g_nesc++] = t;
}

/* Read the sph file once, classifying each line into a dot (resid -999), a
   keep line (everything else except -888), or dropped (-888, or malformed
   coordinates) - the exact filter _conn_classify_sph applies, in file order
   (later index semantics for lobe membership depend on this order). */
static void read_sph(const char *path, double ux, double uy, double uz, double ox, double oy, double oz) {
    hio_reader rd; hio_rec a;
    if (!hio_open(&rd, path)) { fprintf(stderr, "cannot read %s\n", path); exit(1); }
    while (hio_next(&rd, &a)) {
        if (a.resseq_ok && a.resseq == -888) continue;
        if (!a.xyz_ok) continue;
        double t = (a.x - ox) * ux + (a.y - oy) * uy + (a.z - oz) * uz;
        if (a.resseq_ok && a.resseq == -999) {
            push_dot(t, a.x, a.y, a.z, a.line);
            if (a.marked) mark_last(0, g_ndots - 1);
            continue;
        }
        push_keep(a.line);
        if (a.beta_ok) {
            if (a.beta > 900.0) push_esc(t);
            else if (a.beta > 0.005) push_cen(t, a.x, a.y, a.z, a.beta);
        }
        if (a.marked) mark_last(1, g_nkeep - 1);
    }
    hio_close(&rd);
}

static int cmp_cen(const void *a, const void *b) {
    double d = ((const Cen *)a)->t - ((const Cen *)b)->t;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}
static int cmp_double(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b;
    return d < 0 ? -1 : d > 0 ? 1 : 0;
}

/* _conn_centreline_at: binary search + linear interpolation on g_cen. */
static void centreline_at(double t, double *cx, double *cy, double *cz, double *wall) {
    const Cen *f = &g_cen[0], *l = &g_cen[g_ncen - 1];
    if (t <= f->t) { *cx = f->x; *cy = f->y; *cz = f->z; *wall = f->r; return; }
    if (t >= l->t) { *cx = l->x; *cy = l->y; *cz = l->z; *wall = l->r; return; }
    int lo = 0, hi = g_ncen - 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        if (g_cen[mid].t <= t) lo = mid; else hi = mid;
    }
    const Cen *a = &g_cen[lo], *bb = &g_cen[hi];
    double d = bb->t - a->t;
    if (d <= 1e-12) { *cx = bb->x; *cy = bb->y; *cz = bb->z; *wall = bb->r; return; }
    double w = (t - a->t) / d;
    *cx = a->x + (bb->x - a->x) * w;
    *cy = a->y + (bb->y - a->y) * w;
    *cz = a->z + (bb->z - a->z) * w;
    *wall = a->r + (bb->r - a->r) * w;
}

/* _conn_axis_basis */
static void axis_basis(double ux, double uy, double uz, double *f1x, double *f1y, double *f1z,
                        double *f2x, double *f2y, double *f2z) {
    double ax, ay, az;
    if (fabs(ux) < 0.9) { ax = 1; ay = 0; az = 0; } else { ax = 0; ay = 1; az = 0; }
    double e1x = ay*uz - az*uy, e1y = az*ux - ax*uz, e1z = ax*uy - ay*ux;
    double n = sqrt(e1x*e1x + e1y*e1y + e1z*e1z);
    if (n < 1e-9) { *f1x=1;*f1y=0;*f1z=0; *f2x=0;*f2y=1;*f2z=0; return; }
    e1x/=n; e1y/=n; e1z/=n;
    *f1x=e1x; *f1y=e1y; *f1z=e1z;
    *f2x = uy*e1z - uz*e1y; *f2y = uz*e1x - ux*e1z; *f2z = ux*e1y - uy*e1x;
}

/* -------------------------------------------------------------- classify */

typedef struct { int idx; double t, az, rr, wall; } LatPt;
typedef struct { int is_pore; LatPt lat; } Cls;   /* per-dot outcome, dot index order preserved */

typedef struct { double z, a; int n; double ef, ka, kb; int has_ka, has_kb; int *members; int nmem; } Lobe;

static int cmp_lobe_n_desc(const void *a, const void *b) {
    return ((const Lobe *)b)->n - ((const Lobe *)a)->n;
}

static void do_classify(double ux, double uy, double uz, double margin,
                        double f1x, double f1y, double f1z, double f2x, double f2y, double f2z) {
    int *pore_dot = xmalloc(g_ndots * sizeof(int));      /* 1 = pore, 0 = lateral */
    LatPt *lat = xmalloc(g_ndots * sizeof(LatPt));        /* lateral dots, in dot order */
    int nlat = 0;
    for (int i = 0; i < g_ndots; i++) {
        double cx, cy, cz, wall;
        centreline_at(g_dots[i].t, &cx, &cy, &cz, &wall);
        double dx = g_dots[i].x - cx, dy = g_dots[i].y - cy, dz = g_dots[i].z - cz;
        double axc = dx*ux + dy*uy + dz*uz;
        dx -= axc*ux; dy -= axc*uy; dz -= axc*uz;
        double rr = sqrt(dx*dx + dy*dy + dz*dz);
        if (rr <= wall + margin) { pore_dot[i] = 1; continue; }
        pore_dot[i] = 0;
        lat[nlat].idx = i;
        lat[nlat].t = g_dots[i].t;
        lat[nlat].az = atan2(dx*f2x + dy*f2y + dz*f2z, dx*f1x + dy*f1y + dz*f1z);
        lat[nlat].rr = rr; lat[nlat].wall = wall;
        nlat++;
    }

    /* escaped ranges: sort, merge within a 2.0 gap */
    qsort(g_esct, g_nesc, sizeof(double), cmp_double);
    double *rlo = xmalloc((g_nesc ? g_nesc : 1) * sizeof(double));
    double *rhi = xmalloc((g_nesc ? g_nesc : 1) * sizeof(double));
    int nranges = 0;
    if (g_nesc) {
        double lo = g_esct[0], hi = lo;
        for (int i = 1; i < g_nesc; i++) {
            if (g_esct[i] - hi > 2.0) { rlo[nranges] = lo; rhi[nranges] = hi; nranges++; lo = g_esct[i]; }
            hi = g_esct[i];
        }
        rlo[nranges] = lo; rhi[nranges] = hi; nranges++;
    }

    /* lobe clustering: grid cells (zcell=3.0, nth=18), flood fill */
    double zcell = 3.0; int nth = 18;
    Lobe *lobes = NULL; int nlobes = 0, lobecap = 0;
    if (nlat >= 20) {
        int zi_min = INT_MAX, zi_max = INT_MIN;
        int *zi = xmalloc(nlat * sizeof(int)), *ti = xmalloc(nlat * sizeof(int));
        for (int i = 0; i < nlat; i++) {
            zi[i] = (int)floor(lat[i].t / zcell);
            int tt = (int)floor((lat[i].az + M_PI) / (2*M_PI) * nth) % nth;
            if (tt < 0) tt += nth;
            ti[i] = tt;
            if (zi[i] < zi_min) zi_min = zi[i];
            if (zi[i] > zi_max) zi_max = zi[i];
        }
        int nz = zi_max - zi_min + 1;
        /* cell -> list of lateral-dot indices (into lat[]) */
        int **cell = xmalloc((size_t)nz * nth * sizeof(int *));
        int *cellcap = xmalloc((size_t)nz * nth * sizeof(int));
        int *celln = xmalloc((size_t)nz * nth * sizeof(int));
        for (int c = 0; c < nz * nth; c++) { cell[c] = NULL; cellcap[c] = 0; celln[c] = 0; }
        for (int i = 0; i < nlat; i++) {
            int c = (zi[i] - zi_min) * nth + ti[i];
            if (celln[c] >= cellcap[c]) { cellcap[c] = cellcap[c] ? cellcap[c]*2 : 8; cell[c] = xrealloc(cell[c], cellcap[c]*sizeof(int)); }
            cell[c][celln[c]++] = i;
        }
        char *seen = xmalloc((size_t)nz * nth); memset(seen, 0, (size_t)nz * nth);
        int *queue = xmalloc((size_t)nz * nth * sizeof(int));
        int *members = xmalloc((size_t)nz * nth * sizeof(int));   /* cell ids in this component */
        /* A near-empty cell does NOT conduct. Measured on a real frame: the
           widest lobe held 1500 dots over 45 cells at a median of 19 per cell,
           but 10 of them held 1-3 - and dropping just those split it in two.
           Those are stray dots, not a passage, and under plain 8-connectivity
           they welded two openings 160 degrees apart into one coloured region.
           Threshold scales with the cloud's own density (dot density is a user
           setting): a fifth of the median occupancy, never below 2.
           Mirrors _conn_frame_lobes in the Tcl port - both must agree. */
        int cond = 2;
        {
            int nocc = 0;
            for (int c = 0; c < nz * nth; c++) if (celln[c] > 0) nocc++;
            if (nocc > 0) {
                int *occs = xmalloc((size_t)nocc * sizeof(int));
                int k = 0;
                for (int c = 0; c < nz * nth; c++) if (celln[c] > 0) occs[k++] = celln[c];
                for (int a = 1; a < nocc; a++) {      /* insertion sort: nocc is small */
                    int v = occs[a], b = a - 1;
                    while (b >= 0 && occs[b] > v) { occs[b+1] = occs[b]; b--; }
                    occs[b+1] = v;
                }
                int med = occs[nocc/2];
                int t = (int)ceil(med * 0.2);
                if (t > cond) cond = t;
                free(occs);
            }
        }
        for (int c0 = 0; c0 < nz * nth; c0++) {
            if (celln[c0] < cond || seen[c0]) continue;
            seen[c0] = 1;
            int qh = 0, qt = 0, mh = 0;
            queue[qt++] = c0;
            while (qh < qt) {
                int c = queue[qh++];
                members[mh++] = c;
                int czi = c / nth, cti = c % nth;
                for (int dz = -1; dz <= 1; dz++) for (int dt = -1; dt <= 1; dt++) {
                    int nzi = czi + dz; if (nzi < 0 || nzi >= nz) continue;
                    int nti = ((cti + dt) % nth + nth) % nth;
                    int nb = nzi * nth + nti;
                    if (celln[nb] >= cond && !seen[nb]) { seen[nb] = 1; queue[qt++] = nb; }
                }
            }
            /* Sparse neighbours join this component rather than being dropped,
               so no dot leaves the picture - they simply cannot BRIDGE two. */
            int mh0 = mh;
            for (int m = 0; m < mh0; m++) {
                int czi = members[m] / nth, cti = members[m] % nth;
                for (int dz = -1; dz <= 1; dz++) for (int dt = -1; dt <= 1; dt++) {
                    int nzi = czi + dz; if (nzi < 0 || nzi >= nz) continue;
                    int nti = ((cti + dt) % nth + nth) % nth;
                    int nb = nzi * nth + nti;
                    if (celln[nb] > 0 && celln[nb] < cond && !seen[nb]) {
                        seen[nb] = 1; members[mh++] = nb;
                    }
                }
            }
            int total = 0;
            for (int m = 0; m < mh; m++) total += celln[members[m]];
            if (total < nlat * 0.02) continue;
            int *idx = xmalloc(total * sizeof(int)); int ni = 0;
            for (int m = 0; m < mh; m++) for (int k = 0; k < celln[members[m]]; k++) idx[ni++] = cell[members[m]][k];
            double zs = 0, sa = 0, ca = 0; int nesc = 0;
            double exit_rr = -1, exit_wall = 0, far_rr = -1, far_wall = 0; int have_exit = 0, have_far = 0;
            for (int q = 0; q < ni; q++) {
                LatPt *p = &lat[idx[q]];
                zs += p->t; sa += sin(p->az); ca += cos(p->az);
                int esc = 0;
                for (int r = 0; r < nranges; r++) if (p->t >= rlo[r]-1.5 && p->t <= rhi[r]+1.5) { esc = 1; break; }
                if (esc) nesc++;
                if (!have_exit || p->rr < exit_rr) { exit_rr = p->rr; exit_wall = p->wall; have_exit = 1; }
                if (!have_far  || p->rr > far_rr)  { far_rr  = p->rr; far_wall  = p->wall; have_far  = 1; }
            }
            if (nlobes >= lobecap) { lobecap = lobecap ? lobecap*2 : 16; lobes = xrealloc(lobes, lobecap*sizeof(Lobe)); }
            Lobe *lb = &lobes[nlobes++];
            lb->z = zs / ni; lb->a = atan2(sa, ca); lb->n = ni;
            lb->ef = (double)nesc / ni;
            lb->has_ka = have_exit; lb->ka = have_exit ? (exit_rr - exit_wall < 0 ? 0.0 : exit_rr - exit_wall) : 0.0;
            lb->has_kb = have_far;
            if (have_far) { double v = far_rr - far_wall - margin; lb->kb = v < 0 ? 0.0 : v; } else lb->kb = 0.0;
            lb->members = idx; lb->nmem = ni;
        }
        qsort(lobes, nlobes, sizeof(Lobe), cmp_lobe_n_desc);
        free(zi); free(ti); free(seen); free(queue); free(members);
        for (int c = 0; c < nz*nth; c++) free(cell[c]);
        free(cell); free(cellcap); free(celln);
    }

    /* ---- output ---- */
    int npore = 0; for (int i = 0; i < g_ndots; i++) npore += pore_dot[i];
    printf("PORE %d\n", npore);
    for (int i = 0; i < g_ndots; i++) if (pore_dot[i]) printf("%s\n", g_dots[i].line);
    printf("KEEP %d\n", g_nkeep);
    for (int i = 0; i < g_nkeep; i++) printf("%s\n", g_keep[i]);
    printf("LATERAL %d\n", nlat);
    for (int i = 0; i < nlat; i++) printf("%s\n", g_dots[lat[i].idx].line);
    printf("ESCRANGE %d\n", nranges);
    for (int i = 0; i < nranges; i++) printf("%.6f %.6f\n", rlo[i], rhi[i]);
    printf("LOBE %d\n", nlobes);
    for (int i = 0; i < nlobes; i++) {
        Lobe *lb = &lobes[i];
        printf("%.6f %.6f %d %.6f %s %s", lb->z, lb->a, lb->n, lb->ef,
               lb->has_ka ? "" : "-", lb->has_kb ? "" : "-");
        if (lb->has_ka) printf("%.6f", lb->ka); else printf("-");
        printf(" ");
        if (lb->has_kb) printf("%.6f", lb->kb); else printf("-");
        for (int m = 0; m < lb->nmem; m++) printf(" %d", lb->members[m]);
        printf("\n");
    }
    int nmarked = 0;
    for (int i = 0; i < g_ndots; i++) nmarked += g_dotmark[i];
    for (int i = 0; i < g_nkeep; i++) nmarked += g_keepmark[i];
    printf("MARKED %d\n", nmarked);
    for (int i = 0; i < g_ndots; i++) if (g_dotmark[i]) printf("%s\n", g_dots[i].line);
    for (int i = 0; i < g_nkeep; i++) if (g_keepmark[i]) printf("%s\n", g_keep[i]);
    fflush(stdout);

    free(pore_dot); free(lat); free(rlo); free(rhi);
    for (int i = 0; i < nlobes; i++) free(lobes[i].members);
    free(lobes);
    for (int i = 0; i < g_ndots; i++) free(g_dots[i].line);
    free(g_dots);
    for (int i = 0; i < g_nkeep; i++) free(g_keep[i]);
    free(g_keep);
    free(g_cen); free(g_esct); free(g_dotmark); free(g_keepmark);
}

/* ------------------------------------------------------------------ split */

typedef struct { double x, y, z; int region; } LDot;

static unsigned cellkey_hash(long a, long b, long c, unsigned mod) {
    unsigned long h = (unsigned long)(a * 73856093L) ^ (unsigned long)(b * 19349663L) ^ (unsigned long)(c * 83492791L);
    return (unsigned)(h % mod);
}

typedef struct HCell { long a, b, c; LDot *dots; int n, cap; int solo; int has_solo; struct HCell *next; } HCell;

typedef struct {
    HCell **buckets; unsigned nb;
} Grid;

static Grid grid_new(unsigned nb) {
    Grid g; g.nb = nb; g.buckets = xmalloc(nb * sizeof(HCell *));
    for (unsigned i = 0; i < nb; i++) g.buckets[i] = NULL;
    return g;
}
static HCell *grid_find(Grid *g, long a, long b, long c, int create) {
    unsigned h = cellkey_hash(a, b, c, g->nb);
    for (HCell *c2 = g->buckets[h]; c2; c2 = c2->next)
        if (c2->a == a && c2->b == b && c2->c == c) return c2;
    if (!create) return NULL;
    HCell *c2 = xmalloc(sizeof(HCell));
    c2->a = a; c2->b = b; c2->c = c; c2->dots = NULL; c2->n = 0; c2->cap = 0;
    c2->has_solo = 0; c2->solo = -1;
    c2->next = g->buckets[h]; g->buckets[h] = c2;
    return c2;
}
static void grid_add(Grid *g, long a, long b, long c, double x, double y, double z, int region) {
    HCell *cell = grid_find(g, a, b, c, 1);
    if (cell->n >= cell->cap) { cell->cap = cell->cap ? cell->cap * 2 : 4; cell->dots = xrealloc(cell->dots, cell->cap * sizeof(LDot)); }
    cell->dots[cell->n].x = x; cell->dots[cell->n].y = y; cell->dots[cell->n].z = z; cell->dots[cell->n].region = region;
    cell->n++;
}

#define REGION_MAX 64
static char *g_region_names[REGION_MAX];
static FILE *g_region_out[REGION_MAX];
static long g_region_count[REGION_MAX];
static int g_nregions = 0;

static void load_labels(Grid *g, const char *path, int rid, double H) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot read %s\n", path); exit(1); }
    char line[512];
    while (fgets(line, sizeof line, f)) {
        size_t len = strlen(line);
        while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
        if (!is_atom_line(line)) continue;
        char xf[16], yf[16], zf[16];
        field(line, len, 30, 37, xf, sizeof xf);
        field(line, len, 38, 45, yf, sizeof yf);
        field(line, len, 46, 53, zf, sizeof zf);
        double x, y, z;
        if (!is_double(xf, &x) || !is_double(yf, &y) || !is_double(zf, &z)) continue;
        long a = (long)floor(x / H), b = (long)floor(y / H), c = (long)floor(z / H);
        grid_add(g, a, b, c, x, y, z, rid);
    }
    fclose(f);
}

static void resolve_solo(Grid *g) {
    for (unsigned i = 0; i < g->nb; i++)
        for (HCell *c = g->buckets[i]; c; c = c->next) {
            int nm = c->dots[0].region, mixed = 0;
            for (int k = 1; k < c->n; k++) if (c->dots[k].region != nm) { mixed = 1; break; }
            if (!mixed) { c->has_solo = 1; c->solo = nm; }
        }
}

/* Nearest region within a growing cubic box, memoized per resolved cell -
   same semantics as _split_conn_mesh_by_region: once a cell (solo or
   resolved-by-search) has an answer, every triangle whose centroid falls in
   it reuses that answer without a fresh search. */
static int classify_centroid(Grid *g, double H, double cx, double cy, double cz) {
    long bi = (long)floor(cx / H), bj = (long)floor(cy / H), bk = (long)floor(cz / H);
    HCell *here = grid_find(g, bi, bj, bk, 0);
    if (here && here->has_solo) return here->solo;
    if (!here) here = grid_find(g, bi, bj, bk, 1);
    if (here->has_solo) return here->solo;   /* resolved by an earlier triangle in this cell */

    int best = -1; double bd = 1e30;
    /* rad grows the box; once ANY dot has been found in a fully-scanned box
       the answer is final, so only the NEW shell at each radius needs
       scanning once the inner box is confirmed empty. */
    for (int rad = 1; rad <= 4 && best < 0; rad++) {
        for (long di = -rad; di <= rad; di++)
        for (long dj = -rad; dj <= rad; dj++)
        for (long dk = -rad; dk <= rad; dk++) {
            if (rad > 1 && labs(di) < rad && labs(dj) < rad && labs(dk) < rad) continue; /* already scanned */
            HCell *c = grid_find(g, bi+di, bj+dj, bk+dk, 0);
            if (!c) continue;
            for (int k = 0; k < c->n; k++) {
                double dx = cx - c->dots[k].x, dy = cy - c->dots[k].y, dz = cz - c->dots[k].z;
                double e = dx*dx + dy*dy + dz*dz;
                if (e < bd) { bd = e; best = c->dots[k].region; }
            }
        }
    }
    here->has_solo = 1; here->solo = best;   /* memoize, matching the Tcl solo() cache */
    return best;
}

static void do_split(const char *union_plot, char **names, char **label_paths, char **out_paths, int nreg) {
    double H = 2.0;
    Grid g = grid_new(65536);
    for (int i = 0; i < nreg; i++) {
        g_region_names[i] = names[i];
        load_labels(&g, label_paths[i], i, H);
    }
    g_nregions = nreg;
    resolve_solo(&g);

    /* "-" registers a region's dots in the classification grid (so nearby
       triangles resolve correctly) without writing its triangles anywhere -
       for a region Tcl already has a valid cached surface for and is not
       rebuilding this call. */
    for (int i = 0; i < nreg; i++) {
        if (!strcmp(out_paths[i], "-")) { g_region_out[i] = NULL; g_region_count[i] = 0; continue; }
        g_region_out[i] = fopen(out_paths[i], "w");
        if (!g_region_out[i]) { fprintf(stderr, "cannot write %s\n", out_paths[i]); exit(1); }
        fprintf(g_region_out[i], "draw delete all\n");
        g_region_count[i] = 0;
    }

    FILE *f = fopen(union_plot, "r");
    if (!f) { fprintf(stderr, "cannot read %s\n", union_plot); exit(1); }
    char *line = NULL; size_t cap = 0; ssize_t n;
    while ((n = getline(&line, &cap, f)) >= 0) {
        while (n > 0 && (line[n-1] == '\n' || line[n-1] == '\r')) line[--n] = 0;
        if (strncmp(line, "draw trinorm ", 13)) continue;
        double v[9]; int nv = 0;
        char *p = line;
        for (int vtx = 0; vtx < 3 && nv < 9; vtx++) {
            p = strchr(p, '{'); if (!p) break; p++;
            char *close = strchr(p, '}'); if (!close) break;
            char buf[128]; size_t blen = (size_t)(close - p); if (blen >= sizeof buf) blen = sizeof buf - 1;
            memcpy(buf, p, blen); buf[blen] = 0;
            char *q = buf;
            v[nv++] = strtod(q, &q); v[nv++] = strtod(q, &q); v[nv++] = strtod(q, &q);
            p = close + 1;
        }
        if (nv < 9) continue;
        double cx = (v[0]+v[3]+v[6]) / 3.0, cy = (v[1]+v[4]+v[7]) / 3.0, cz = (v[2]+v[5]+v[8]) / 3.0;
        int rid = classify_centroid(&g, H, cx, cy, cz);
        if (rid < 0 || rid >= nreg) continue;
        g_region_count[rid]++;
        if (!g_region_out[rid]) continue;
        fprintf(g_region_out[rid], "%s\n", line);
    }
    free(line);
    fclose(f);
    long total = 0;
    for (int i = 0; i < nreg; i++) { if (g_region_out[i]) fclose(g_region_out[i]); total += g_region_count[i]; }
    printf("%ld\n", total);

    for (unsigned i = 0; i < g.nb; i++) {
        HCell *c = g.buckets[i];
        while (c) { HCell *next = c->next; free(c->dots); free(c); c = next; }
    }
    free(g.buckets);
}

/* ------------------------------------------------------------------ main */

static void usage(const char *a0) {
    fprintf(stderr,
        "usage: %s classify SPH CX CY CZ VX VY VZ MARGIN [F1X F1Y F1Z F2X F2Y F2Z]\n"
        "       %s split UNION_PLOT --region NAME LABELS_SPH OUT_PLOT [...]\n",
        a0, a0);
}

#ifdef VMDPATHFINDER_MULTICALL
int conn_lobes_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
    if (argc < 2) { usage(argv[0]); return 2; }
    if (!strcmp(argv[1], "classify")) {
        if (argc < 10) { usage(argv[0]); return 2; }
        const char *sph = argv[2];
        double cx = atof(argv[3]), cy = atof(argv[4]), cz = atof(argv[5]);
        double vx = atof(argv[6]), vy = atof(argv[7]), vz = atof(argv[8]);
        double margin = atof(argv[9]);
        double ulen = sqrt(vx*vx + vy*vy + vz*vz);
        if (ulen <= 1e-9) { fprintf(stderr, "CVECT is a zero vector\n"); return 1; }
        double ux = vx/ulen, uy = vy/ulen, uz = vz/ulen;
        double f1x, f1y, f1z, f2x, f2y, f2z;
        if (argc >= 16) {
            f1x=atof(argv[10]); f1y=atof(argv[11]); f1z=atof(argv[12]);
            f2x=atof(argv[13]); f2y=atof(argv[14]); f2z=atof(argv[15]);
        } else {
            axis_basis(ux, uy, uz, &f1x, &f1y, &f1z, &f2x, &f2y, &f2z);
        }
        read_sph(sph, ux, uy, uz, cx, cy, cz);
        if (g_ncen < 2 || g_ndots == 0) { printf("PORE 0\nKEEP 0\nLATERAL 0\nESCRANGE 0\nLOBE 0\n"); return 0; }
        qsort(g_cen, g_ncen, sizeof(Cen), cmp_cen);
        do_classify(ux, uy, uz, margin, f1x, f1y, f1z, f2x, f2y, f2z);
        return 0;
    }
    if (!strcmp(argv[1], "split")) {
        if (argc < 4) { usage(argv[0]); return 2; }
        const char *union_plot = argv[2];
        char *names[REGION_MAX], *labels[REGION_MAX], *outs[REGION_MAX];
        int nreg = 0;
        int i = 3;
        while (i < argc) {
            if (strcmp(argv[i], "--region") || i + 3 >= argc) { usage(argv[0]); return 2; }
            if (nreg >= REGION_MAX) { fprintf(stderr, "too many regions\n"); return 2; }
            names[nreg] = argv[i+1]; labels[nreg] = argv[i+2]; outs[nreg] = argv[i+3];
            nreg++; i += 4;
        }
        if (!nreg) { usage(argv[0]); return 2; }
        do_split(union_plot, names, labels, outs, nreg);
        return 0;
    }
    usage(argv[0]);
    return 2;
}
