/* MOLE 2's tunnel stage: cavity openings, Dijkstra, profile, filters.
 *
 * Ported from WebChemistry.Tunnels.Core (MIT). Checkpoint for this stage is
 * 3 tunnels on 1tqn with test.xml's origin.
 *
 * Note the SurfaceCavity contributes nothing here: Cavity.CreateSurface leaves
 * Openings empty, and GetTunnels returns early for a source with no openings, so
 * it only matters once custom exits exist. Only the origin's own cavity is
 * searched.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_tunnel.h"
#include "../xalloc.h"

/* ---- openings ---------------------------------------------------------- */

/* CavityOpening.Common: two boundary tetrahedra are joined when they share at
   least TWO atoms - not when they are face-adjacent. A weaker relation, so the
   components are coarser than the cavity graph's own. */
static int shares_two_atoms(const mole_complex *M, int a, int b)
{
    int i, j, cv = 0;
    for (i = 0; i < 4; i++) {
        for (j = 0; j < 4; j++)
            if (M->tv[4*a+i] == M->tv[4*b+j]) { cv++; break; }
        if (cv > 1) return 1;
    }
    return 0;
}

/* CavityOpening.Cover: greedily take the widest tetrahedron as a pivot, discard
   everything within coverRadius of its CENTRE, recurse on what is left. This is
   what SurfaceCoverRadius controls, and it is why one boundary component can
   yield several openings. */
static void cover(const mole_complex *M, int *comp, int n, double radius,
                  int *pivots, int *npiv)
{
    int i, best = -1, m = 0;
    int *rest;
    if (n <= 0) return;
    for (i = 0; i < n; i++)
        if (best < 0 || M->maxclear[comp[i]] > M->maxclear[best]) best = comp[i];
    pivots[(*npiv)++] = best;

    rest = malloc((size_t)n * sizeof(int));
    if (!rest) return;
    for (i = 0; i < n; i++) {
        double dx = M->center[3*comp[i]]   - M->center[3*best];
        double dy = M->center[3*comp[i]+1] - M->center[3*best+1];
        double dz = M->center[3*comp[i]+2] - M->center[3*best+2];
        if (sqrt(dx*dx + dy*dy + dz*dz) > radius) rest[m++] = comp[i];
    }
    if (m > 0) cover(M, rest, m, radius, pivots, npiv);
    free(rest);
}

int mole_openings(const mole_complex *M, int cavity_comp, double cover_radius,
                  int **out)
{
    int i, j, nb = 0, npiv = 0, *bnd, *lab, *pivots, ncomp = 0;
    for (i = 0; i < M->nt; i++)
        if (M->alive[i] && M->comp[i] == cavity_comp && M->boundary[i]) nb++;
    if (!nb) { *out = NULL; return 0; }

    bnd = malloc((size_t)nb * sizeof(int));
    lab = malloc((size_t)nb * sizeof(int));
    pivots = malloc((size_t)nb * sizeof(int));
    if (!bnd || !lab || !pivots) { free(bnd); free(lab); free(pivots); *out=NULL; return 0; }
    nb = 0;
    for (i = 0; i < M->nt; i++)
        if (M->alive[i] && M->comp[i] == cavity_comp && M->boundary[i]) bnd[nb++] = i;
    /* Cavity.cs:235: boundaryTetra = graph.Vertices.Where(v => v.IsBoundary),
       so candidates arrive in cavityGraph.Vertices order, not by index. That
       order resolves MaxBy's tie (Create takes [0] of all maxima), Cover's
       recursion and the component ids. */
    for (i = 1; i < nb; i++) {
        int v = bnd[i], j = i - 1;
        while (j >= 0 && M->vorder[bnd[j]] > M->vorder[v]) { bnd[j+1] = bnd[j]; j--; }
        bnd[j+1] = v;
    }
    for (i = 0; i < nb; i++) lab[i] = -1;

    /* Connected components under the shares-two-atoms relation. */
    for (i = 0; i < nb; i++) {
        int *stack, sp = 0;
        if (lab[i] >= 0) continue;
        stack = malloc((size_t)nb * sizeof(int));
        if (!stack) break;
        lab[i] = ncomp; stack[sp++] = i;
        while (sp > 0) {
            int c = stack[--sp];
            for (j = 0; j < nb; j++)
                if (lab[j] < 0 && shares_two_atoms(M, bnd[c], bnd[j])) {
                    lab[j] = ncomp; stack[sp++] = j;
                }
        }
        free(stack);
        ncomp++;
    }

    for (i = 0; i < ncomp; i++) {
        int *members = malloc((size_t)nb * sizeof(int)), nm = 0;
        if (!members) break;
        for (j = 0; j < nb; j++) if (lab[j] == i) members[nm++] = bnd[j];
        cover(M, members, nm, cover_radius, pivots, &npiv);
        free(members);
    }
    free(bnd); free(lab);
    *out = pivots;
    return npiv;
}

/* TunnelOriginCollection.FromCavity: the computed start points.
 *
 * Both method papers define the start point by the INTEGER depth h. The code
 * uses real-valued DepthLength as primary, with the integer version only as a
 * fallback, and adds two conditions that appear in neither paper: a floor at a
 * quarter of the cavity's maximum DepthLength, and a minimum degree of 3.
 * Candidates are local maxima of DepthLength, taken deepest-first and greedily
 * spread by AutoOriginCoverRadius up to MaxAutoOriginsPerCavity.
 */
int mole_auto_origins(const mole_complex *M, int cavity_comp,
                      double cover_radius, int max_origins, int *out)
{
    int i, k, nc = 0, n = 0, *cand;
    double mx = 0.0, floor_;

    /* `out` is a fixed MOLE_MAX_ORIGINS stack array in the caller. Clamp here
       too: the caller already clamps, but this function is the one that
       writes, and the unconditional out[n++] below would overflow even for
       max_origins == 0. */
    if (max_origins > MOLE_MAX_ORIGINS) max_origins = MOLE_MAX_ORIGINS;
    if (max_origins < 1) return 0;

    for (i = 0; i < M->nt; i++)
        if (M->alive[i] && M->comp[i] == cavity_comp
            && M->depthlen[i] < 1e299 && M->depthlen[i] > mx) mx = M->depthlen[i];
    floor_ = mx / 4.0;

    cand = malloc((size_t)M->nt * sizeof(int));
    if (!cand) return 0;
    for (i = 0; i < M->nt; i++) {
        int deg = 0, is_max = 1;
        if (!M->alive[i] || M->comp[i] != cavity_comp) continue;
        if (!(M->depthlen[i] > floor_)) continue;
        for (k = 0; k < 4; k++) {
            int v = M->tn[4*i+k];
            if (v < 0 || !M->alive[v] || M->comp[v] != cavity_comp) continue;
            deg++;
            if (M->depthlen[v] >= M->depthlen[i]) is_max = 0;
        }
        if (is_max && deg >= 3) cand[nc++] = i;
    }
    /* descending DepthLength; insertion sort, the candidate list is short */
    for (i = 1; i < nc; i++) {
        int v = cand[i], j = i - 1;
        while (j >= 0 && M->depthlen[cand[j]] < M->depthlen[v]) { cand[j+1] = cand[j]; j--; }
        cand[j+1] = v;
    }
    /* FALLBACK, and it is the common case: on 1AKD the primary list is empty in
       four of five cavities because the deepest vertex has degree 1. MOLE then
       uses integer-Depth local maxima, with NO floor and NO degree condition,
       keeps only those attaining the maximum depth, and spreads them by
       AutoOriginCoverRadius.

       Note the spreading test compares each candidate against all EARLIER
       CANDIDATES, not against the ones actually added - which is not the same
       greedy, and is written the way they write it. */
    if (nc == 0) {
        int mxd = -1;
        for (i = 0; i < M->nt; i++) {
            int is_max = 1, has_nb = 0;
            if (!M->alive[i] || M->comp[i] != cavity_comp) continue;
            for (k = 0; k < 4; k++) {
                int v = M->tn[4*i+k];
                if (v < 0 || !M->alive[v] || M->comp[v] != cavity_comp) continue;
                has_nb = 1;
                if (M->depth[v] >= M->depth[i]) is_max = 0;
            }
            if (has_nb && is_max) { cand[nc++] = i; if (M->depth[i] > mxd) mxd = M->depth[i]; }
        }
        if (nc == 0) { free(cand); return 0; }
        { int m2 = 0;
          for (i = 0; i < nc; i++) if (M->depth[cand[i]] == mxd) cand[m2++] = cand[i];
          nc = m2; }
        for (i = 0; i < nc && n < max_origins; i++) {
            int add = 1, j;
            for (j = i - 1; j >= 0; j--) {
                double dx = M->center[3*cand[i]]   - M->center[3*cand[j]];
                double dy = M->center[3*cand[i]+1] - M->center[3*cand[j]+1];
                double dz = M->center[3*cand[i]+2] - M->center[3*cand[j]+2];
                if (sqrt(dx*dx+dy*dy+dz*dz) < cover_radius) { add = 0; break; }
            }
            if (add) out[n++] = cand[i];
        }
        free(cand);
        return n;
    }

    out[n++] = cand[0];   /* safe: max_origins >= 1 checked on entry */
    for (;;) {
        int added = 0;
        for (i = 0; i < nc && n < max_origins; i++) {
            int ok = 1, j;
            for (j = 0; j < n; j++) {
                double dx = M->center[3*cand[i]]   - M->center[3*out[j]];
                double dy = M->center[3*cand[i]+1] - M->center[3*out[j]+1];
                double dz = M->center[3*cand[i]+2] - M->center[3*out[j]+2];
                if (sqrt(dx*dx+dy*dy+dz*dz) <= cover_radius) { ok = 0; break; }
            }
            if (ok) { out[n++] = cand[i]; added = 1; break; }
        }
        if (!added || n >= max_origins) break;
    }
    free(cand);
    return n;
}

/* ---- Dijkstra ---------------------------------------------------------- */

typedef struct { int v; double d; } heap_item;

static void heap_push(heap_item *h, int *n, int v, double d)
{
    int i = (*n)++;
    h[i].v = v; h[i].d = d;
    while (i > 0) {
        int p = (i-1)/2;
        if (h[p].d <= h[i].d) break;
        { heap_item t = h[p]; h[p] = h[i]; h[i] = t; i = p; }
    }
}

static int heap_pop(heap_item *h, int *n)
{
    int best = h[0].v, i = 0;
    h[0] = h[--(*n)];
    for (;;) {
        int l = 2*i+1, r = l+1, s = i;
        if (l < *n && h[l].d < h[s].d) s = l;
        if (r < *n && h[r].d < h[s].d) s = r;
        if (s == i) break;
        { heap_item t = h[s]; h[s] = h[i]; h[i] = t; i = s; }
    }
    return best;
}

/* Complex.MakeEdgeWeightFunction. */
double mole_edge_cost(const mole_complex *M, int t, int k, mole_weight_fn w)
{
    switch (w) {
    case MOLE_W_LENGTH_AND_RADIUS: return M->eweight[4*t+k];
    case MOLE_W_LENGTH:            return M->elen[4*t+k];
    case MOLE_W_CONSTANT:          return 1.0;
    default:                       return M->evweight[4*t+k];
    }
}

/* mask, when given, REPLACES the alive/comp membership test - that is how the
   SurfaceCavity is searched, since its tetrahedra are not alive by the time the
   tunnel stage runs and belong to no cavity component. NULL keeps the ordinary
   cavity behaviour exactly. */
int mole_dijkstra_mask(const mole_complex *M, int cavity_comp, int src,
                       double *dist, int *prev, mole_weight_fn w,
                       const char *mask)
{
    heap_item *h = malloc((size_t)(M->nt*4+8) * sizeof(*h));
    int nh = 0, i;
    if (!h) return -1;
    for (i = 0; i < M->nt; i++) { dist[i] = 1e300; prev[i] = -1; }
    dist[src] = 0.0;
    heap_push(h, &nh, src, 0.0);
    while (nh > 0) {
        int u = heap_pop(h, &nh), k;
        for (k = 0; k < 4; k++) {
            int v = M->tn[4*u+k];
            double nd;
            if (v < 0) continue;
            if (mask) { if (!mask[v]) continue; }
            else if (!M->alive[v] || M->comp[v] != cavity_comp) continue;
            nd = dist[u] + mole_edge_cost(M, u, k, w);
            if (nd < dist[v]) { dist[v] = nd; prev[v] = u; heap_push(h, &nh, v, nd); }
        }
    }
    free(h);
    return 0;
}

int mole_dijkstra(const mole_complex *M, int cavity_comp, int src,
                  double *dist, int *prev, mole_weight_fn w)
{
    return mole_dijkstra_mask(M, cavity_comp, src, dist, prev, w, NULL);
}

/* ---- natural cubic spline ---------------------------------------------- */

/* Zero second derivative at both ends, which is what MOLE's two-argument
   CubicSplineInterpolation resolves to. Samples are uniform in the parameter t,
   NOT in arc length - the tunnel's own geometry does not reparameterise it. */
void mole_spline_init(mole_spline *s, const double *t, const double *y, int n)
{
    double *a = xa_malloc((size_t)n*sizeof(double));
    double *b = xa_malloc((size_t)n*sizeof(double));
    double *c = xa_malloc((size_t)n*sizeof(double));
    double *d = xa_malloc((size_t)n*sizeof(double));
    int i;
    s->n = n;
    s->t = xa_malloc((size_t)n*sizeof(double));
    s->y = xa_malloc((size_t)n*sizeof(double));
    s->m = xa_malloc((size_t)n*sizeof(double));
    if (!a||!b||!c||!d||!s->t||!s->y||!s->m) return;
    memcpy(s->t, t, (size_t)n*sizeof(double));
    memcpy(s->y, y, (size_t)n*sizeof(double));

    b[0] = 1; c[0] = 0; d[0] = 0;
    for (i = 1; i < n-1; i++) {
        double h0 = t[i]-t[i-1], h1 = t[i+1]-t[i];
        a[i] = h0; b[i] = 2*(h0+h1); c[i] = h1;
        d[i] = 6*((y[i+1]-y[i])/h1 - (y[i]-y[i-1])/h0);
    }
    a[n-1] = 0; b[n-1] = 1; d[n-1] = 0;
    for (i = 1; i < n; i++) {
        double w = a[i]/b[i-1];
        b[i] -= w*c[i-1];
        d[i] -= w*d[i-1];
    }
    s->m[n-1] = d[n-1]/b[n-1];
    for (i = n-2; i >= 0; i--) s->m[i] = (d[i] - c[i]*s->m[i+1])/b[i];
    free(a); free(b); free(c); free(d);
}

double mole_spline_eval(const mole_spline *s, double x)
{
    int lo = 0, hi = s->n-1;
    double h, A, B;
    if (x <= s->t[0]) return s->y[0];
    if (x >= s->t[s->n-1]) return s->y[s->n-1];
    while (hi - lo > 1) { int mid = (lo+hi)/2; if (s->t[mid] > x) hi = mid; else lo = mid; }
    h = s->t[hi]-s->t[lo];
    A = (s->t[hi]-x)/h; B = (x-s->t[lo])/h;
    return A*s->y[lo] + B*s->y[hi]
         + ((A*A*A-A)*s->m[lo] + (B*B*B-B)*s->m[hi]) * (h*h)/6.0;
}

void mole_spline_free(mole_spline *s)
{
    free(s->t); free(s->y); free(s->m);
    memset(s, 0, sizeof(*s));
}


/* Cavity.GetOpening: over the cavity's BOUNDARY facets - which exist only for
   tetrahedra with fewer than four neighbours INSIDE that cavity - take the one
   whose tetrahedron CENTROID (not circumcentre) is nearest the point, and
   reject it past radius. Returns the pivot tetrahedron or -1.

   mask selects the cavity: the surface membership for the SurfaceCavity, NULL
   to use alive/comp as the regular cavities do. */
int mole_cavity_opening(const mole_complex *M, int cavity_comp, const char *mask,
                        const double *point, double radius)
{
    int t, k, best = -1;
    double bd = 1e300;
    for (t = 0; t < M->nt; t++) {
        int deg = 0;
        double dx, dy, dz, d2;
        if (mask) { if (!mask[t]) continue; }
        else if (!M->alive[t] || M->comp[t] != cavity_comp) continue;
        for (k = 0; k < 4; k++) {
            int n = M->tn[4*t+k];
            if (n < 0) continue;
            if (mask) { if (mask[n]) deg++; }
            else if (M->alive[n] && M->comp[n] == cavity_comp) deg++;
        }
        if (deg >= 4) continue;             /* contributes no boundary facet */
        dx = M->center[3*t] - point[0];
        dy = M->center[3*t+1] - point[1];
        dz = M->center[3*t+2] - point[2];
        d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < bd) { bd = d2; best = t; }
    }
    if (best < 0 || bd > radius*radius) return -1;
    return best;
}

/* Cavity.GetTetrahedron: nearest member by CENTROID over the WHOLE cavity -
   unlike GetOpening, which only considers boundary-facet tetrahedra. */
int mole_cavity_tetrahedron(const mole_complex *M, int cavity_comp,
                            const char *mask, const double *point, double radius)
{
    int t, best = -1;
    double bd = 1e300;
    for (t = 0; t < M->nt; t++) {
        double dx, dy, dz, d2;
        if (mask) { if (!mask[t]) continue; }
        else if (!M->alive[t] || M->comp[t] != cavity_comp) continue;
        dx = M->center[3*t] - point[0];
        dy = M->center[3*t+1] - point[1];
        dz = M->center[3*t+2] - point[2];
        d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < bd) { bd = d2; best = t; }
    }
    if (best < 0 || bd > radius*radius) return -1;
    return best;
}

/* ---- profile ----------------------------------------------------------- */

/* Radius at a point: the minimum over the FIVE nearest atoms of
   (distance - vdW), floored at 0.01. Both method papers say ten; the code says
   five, and the code is what produced the reference numbers.
   Brute force over the atom list - correctness first, and the profile is 100
   samples per tunnel. */
/* Uniform grid over the atoms, for the 5-nearest query.
 *
 * radius_at is 36% of runtime as a brute-force scan - 100 samples per tunnel
 * times every atom. The grid must return the SAME five atoms, so the search
 * expands whole cell shells and only stops once the 5th-best distance is
 * provably closer than anything the next shell could hold. Same answer, and
 * verified byte-identical against the four exact structures.
 */
typedef struct {
    double ox, oy, oz, cell;
    int nx, ny, nz;
    int *head, *next;      /* per-cell atom chain */
} atom_grid;

static atom_grid AG;
static const double *AG_xyz;
static int AG_n;

static void grid_build(const double *axyz, int na)
{
    double lo[3], hi[3];
    int i, k;
    if (AG_xyz == axyz && AG_n == na) return;         /* already built */
    free(AG.head); free(AG.next);
    for (k = 0; k < 3; k++) { lo[k] = hi[k] = axyz[k]; }
    for (i = 1; i < na; i++)
        for (k = 0; k < 3; k++) {
            if (axyz[3*i+k] < lo[k]) lo[k] = axyz[3*i+k];
            if (axyz[3*i+k] > hi[k]) hi[k] = axyz[3*i+k];
        }
    AG.cell = 4.0;                                     /* ~2 atoms across */
    AG.ox = lo[0]; AG.oy = lo[1]; AG.oz = lo[2];
    AG.nx = (int)((hi[0]-lo[0])/AG.cell) + 2;
    AG.ny = (int)((hi[1]-lo[1])/AG.cell) + 2;
    AG.nz = (int)((hi[2]-lo[2])/AG.cell) + 2;
    AG.head = xa_malloc((size_t)AG.nx*AG.ny*AG.nz*sizeof(int));
    AG.next = xa_malloc((size_t)na*sizeof(int));
    if (!AG.head || !AG.next) { free(AG.head); free(AG.next); AG.head=NULL; AG.next=NULL; return; }
    for (i = 0; i < AG.nx*AG.ny*AG.nz; i++) AG.head[i] = -1;
    for (i = 0; i < na; i++) {
        int cx = (int)((axyz[3*i]  -AG.ox)/AG.cell);
        int cy = (int)((axyz[3*i+1]-AG.oy)/AG.cell);
        int cz = (int)((axyz[3*i+2]-AG.oz)/AG.cell);
        int c = (cz*AG.ny + cy)*AG.nx + cx;
        AG.next[i] = AG.head[c]; AG.head[c] = i;
    }
    AG_xyz = axyz; AG_n = na;
}

/* sel/nsel, when given, report WHICH five were kept - the lining stage needs the
   atoms themselves, not just the width they imply. */
static double radius_at_impl(const double *p, const double *axyz, const double *arad,
                             int na, int clamp, int *sel, int *nsel)
{
    /* Keep the five nearest by DISTANCE, then take the minimum of
       (distance - vdW) over exactly those five. Selecting directly by
       (distance - vdW) would pick a different, generally larger, set.
       Ranked by SQUARED distance, which is what MOLE's PriorityArray orders on
       (K3DNodes.cs:248, DistanceToSquared). Identical to ranking by the
       distance except when two are within an ULP or so of each other, where
       sqrt can round the comparison the other way - and on a symmetric
       multimer that is the normal case, not a curiosity. bd holds squares
       throughout; the width takes sqrt of the same double, so radii are
       unchanged.

       EXACT TIES break by ATOM INDEX, not by traversal order. The grid path
       below walks atoms in cell-bucket order while the brute-force path walks
       them in index order, so "first encountered wins" made the two C paths
       disagree with EACH OTHER on a tie - and the Tcl port, which scans
       ascending, disagreed with the grid path too. Radii could differ (0.6 vs
       0.2) and lining identities differ even where radii agreed. Exact ties
       are the NORMAL case on a symmetric multimer, not a curiosity. Ranking by
       (squared distance, atom index) is a total order, so all three now give
       the same five atoms in the same sequence regardless of how they were
       reached. This does NOT claim to match real MOLE on a tie: MOLE's own
       choice follows .NET Dictionary iteration order and is unrecoverable
       (see NOTES - the 2OAR tie was traced to exactly that). */
    double bd[5];
    int bi[5], nb = 0, i, j;
    double r = 1e300;

    grid_build(axyz, na);
    if (AG.head) {
        int cx = (int)((p[0]-AG.ox)/AG.cell), cy = (int)((p[1]-AG.oy)/AG.cell);
        int cz = (int)((p[2]-AG.oz)/AG.cell), ring;
        if (cx < 0) cx = 0; if (cx >= AG.nx) cx = AG.nx-1;
        if (cy < 0) cy = 0; if (cy >= AG.ny) cy = AG.ny-1;
        if (cz < 0) cz = 0; if (cz >= AG.nz) cz = AG.nz-1;
        for (ring = 0; ring < AG.nx + AG.ny + AG.nz; ring++) {
            int ax, ay, az;
            /* Stop only when the 5th best is closer than the nearest point the
               next ring could contain - otherwise the answer could change. */
            if (nb == 5 && ring > 0
                && bd[4] <= ((ring - 1) * AG.cell) * ((ring - 1) * AG.cell)) break;
            for (az = cz-ring; az <= cz+ring; az++) {
                if (az < 0 || az >= AG.nz) continue;
                for (ay = cy-ring; ay <= cy+ring; ay++) {
                    if (ay < 0 || ay >= AG.ny) continue;
                    for (ax = cx-ring; ax <= cx+ring; ax++) {
                        int c, q;
                        if (ax < 0 || ax >= AG.nx) continue;
                        /* shell only: skip the interior already visited */
                        if (ring > 0 && abs(ax-cx) < ring && abs(ay-cy) < ring
                            && abs(az-cz) < ring) continue;
                        c = (az*AG.ny + ay)*AG.nx + ax;
                        for (q = AG.head[c]; q >= 0; q = AG.next[q]) {
                            double dx = axyz[3*q]-p[0], dy = axyz[3*q+1]-p[1], dz = axyz[3*q+2]-p[2];
                            double d = dx*dx + dy*dy + dz*dz;
                            if (nb < 5) {
                                bd[nb] = d; bi[nb] = q; nb++;
                                for (j = nb-1; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                                    double td = bd[j]; int ti = bi[j];
                                    bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
                                }
                            } else if (d < bd[4] || (d == bd[4] && q < bi[4])) {
                                bd[4] = d; bi[4] = q;
                                for (j = 4; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                                    double td = bd[j]; int ti = bi[j];
                                    bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
                                }
                            }
                        }
                    }
                }
            }
        }
        for (j = 0; j < nb; j++) { double v = sqrt(bd[j]) - arad[bi[j]]; if (v < r) r = v; }
        if (sel) { for (j = 0; j < nb; j++) sel[j] = bi[j]; *nsel = nb; }
        return (clamp && r < 0.01) ? 0.01 : r;
    }

    for (i = 0; i < na; i++) {
        double dx = axyz[3*i]-p[0], dy = axyz[3*i+1]-p[1], dz = axyz[3*i+2]-p[2];
        double d = dx*dx + dy*dy + dz*dz;
        if (nb < 5) {
            bd[nb] = d; bi[nb] = i; nb++;
            for (j = nb-1; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j] = bd[j-1]; bi[j] = bi[j-1];
                bd[j-1] = td; bi[j-1] = ti;
            }
        } else if (d < bd[4] || (d == bd[4] && i < bi[4])) {
            bd[4] = d; bi[4] = i;
            for (j = 4; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j] = bd[j-1]; bi[j] = bi[j-1];
                bd[j-1] = td; bi[j-1] = ti;
            }
        }
    }
    for (j = 0; j < nb; j++) {
        double v = sqrt(bd[j]) - arad[bi[j]];
        if (v < r) r = v;
    }
    if (sel) { for (j = 0; j < nb; j++) sel[j] = bi[j]; *nsel = nb; }
    return (clamp && r < 0.01) ? 0.01 : r;
}

double mole_radius_at_raw(const double *p, const double *axyz, const double *arad, int na)
{ return radius_at_impl(p, axyz, arad, na, 0, NULL, NULL); }

static double radius_at(const double *p, const double *axyz, const double *arad, int na)
{ return radius_at_impl(p, axyz, arad, na, 1, NULL, NULL); }

/* The five atoms the radius at p was measured against, nearest first. Same
   query, same answer - the lining stage groups them into residues. */
int mole_nearest5(const double *p, const double *axyz, const double *arad,
                  int na, int *sel)
{
    int nb = 0;
    radius_at_impl(p, axyz, arad, na, 1, sel, &nb);
    return nb;
}

/* ---- FreeRadius / BRadius ---------------------------------------------- */

/* Set once per run by the caller; NULL leaves the extra widths equal to the
   plain radius, which is what happens when the atom table carries neither
   column. Globals rather than parameters because mole_profile's signature is
   shared with the Tcl port's trace comparison and every existing caller. */
static const double *PX_BFAC;
static const int    *PX_FREE;

void mole_profile_extras(const double *bfac, const int *freeatom)
{ PX_BFAC = bfac; PX_FREE = freeatom; }

/* radius + the MEAN over the SAME five nearest atoms of the B-factor RMSF,
   sqrt(3B / 8 pi^2). Note MOLE adds the mean to the UNCLAMPED radius and
   clamps the sum, not the other way round. */
static double bradius_at(const double *p, const double *axyz, const double *arad, int na)
{
    double bd[5]; int bi[5], nb = 0, i, j;
    double r = 1e300, s = 0.0;
    if (!PX_BFAC) return radius_at(p, axyz, arad, na);
    for (i = 0; i < na; i++) {
        double dx = axyz[3*i]-p[0], dy = axyz[3*i+1]-p[1], dz = axyz[3*i+2]-p[2];
        double d = dx*dx + dy*dy + dz*dz;
        if (nb < 5) {
            bd[nb] = d; bi[nb] = i; nb++;
            for (j = nb-1; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
            }
        } else if (d < bd[4] || (d == bd[4] && i < bi[4])) {
            bd[4] = d; bi[4] = i;
            for (j = 4; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
            }
        }
    }
    for (j = 0; j < nb; j++) {
        double v = sqrt(bd[j]) - arad[bi[j]];
        if (v < r) r = v;
        s += sqrt(3.0 * PX_BFAC[bi[j]] / (8.0 * 3.14159265358979323846 * 3.14159265358979323846));
    }
    if (nb) s /= nb;
    r += s;
    return r < 0.01 ? 0.01 : r;
}

/* The same query restricted to the atoms FreeRadius is measured against -
   MOLE builds a SECOND k-d tree over backbone and het atoms only, so this is a
   different five, not a filter applied afterwards. */
static double free_radius_at(const double *p, const double *axyz, const double *arad, int na)
{
    double bd[5]; int bi[5], nb = 0, i, j;
    double r = 1e300;
    if (!PX_FREE) return radius_at(p, axyz, arad, na);
    for (i = 0; i < na; i++) {
        double dx, dy, dz, d;
        if (!PX_FREE[i]) continue;
        dx = axyz[3*i]-p[0]; dy = axyz[3*i+1]-p[1]; dz = axyz[3*i+2]-p[2];
        d = dx*dx + dy*dy + dz*dz;
        if (nb < 5) {
            bd[nb] = d; bi[nb] = i; nb++;
            for (j = nb-1; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
            }
        } else if (d < bd[4] || (d == bd[4] && i < bi[4])) {
            bd[4] = d; bi[4] = i;
            for (j = 4; j > 0 && (bd[j] < bd[j-1] || (bd[j] == bd[j-1] && bi[j] < bi[j-1])); j--) {
                double td = bd[j]; int ti = bi[j];
                bd[j]=bd[j-1]; bi[j]=bi[j-1]; bd[j-1]=td; bi[j-1]=ti;
            }
        }
    }
    if (!nb) return radius_at(p, axyz, arad, na);
    for (j = 0; j < nb; j++) {
        double v = sqrt(bd[j]) - arad[bi[j]];
        if (v < r) r = v;
    }
    return r < 0.01 ? 0.01 : r;
}

/* Tetrahedron.ContainsPoint: the barycentric sign test as they write it - four
   determinants with one row replaced by the query point, every sign agreeing
   with the whole, and the sum matching to 1e-6. */
static int contains_point(const mole_complex *M, int t, const double *q)
{
    double m[4][4], d0, s = 0.0, di;
    int i, k;
    for (i = 0; i < 4; i++) {
        int a = M->tv[4*t+i];
        for (k = 0; k < 3; k++) m[i][k] = M->axyz[3*a+k];
        m[i][3] = 1.0;
    }
    d0 = mole_det4(m);
    for (i = 0; i < 4; i++) {
        double save[3];
        for (k = 0; k < 3; k++) { save[k] = m[i][k]; m[i][k] = q[k]; }
        di = mole_det4(m);
        for (k = 0; k < 3; k++) m[i][k] = save[k];
        if (((d0 > 0) - (d0 < 0)) != ((di > 0) - (di < 0))) return 0;
        s += di;
    }
    return fabs(s - d0) > 0.000001 ? 0 : 1;
}

/* Tunnel.CalculateProfile's control path.
 *
 * The spline is NOT fitted to the Dijkstra path. Four steps come first, and
 * omitting them cost 1-4 A of length and one whole tunnel:
 *   - advance past the start while the 5-nearest radius is below the interior
 *     threshold, i.e. while still inside the narrow mouth
 *   - keep path[i] only when its centroid is more than 0.7 A from path[i-1]'s.
 *     The comparison is against the previous PATH element, not the last one
 *     kept, so it filters consecutive originals rather than resampling greedily
 *   - trim from the END while the last tetrahedron neither contains its own
 *     circumcentre nor has it within 3 A of its centroid
 *   - reject the tunnel outright if fewer than 5 control points survive
 */
/* is_path selects TunnelType.Path: CalculateProfile guards the leading skip
   with `if (this.Type == TunnelType.Tunnel)`, so a Path keeps its first
   tetrahedron whatever its radius. */
int mole_control_path_ex(const mole_complex *M, const int *path, int np,
                         const double *axyz, const double *arad, int na,
                         double interior_threshold, int *out, int is_path)
{
    int start, i, n = 0;
    if (is_path) start = 0;
    else {
        for (start = 0; start < np; start++)
            if (mole_radius_at_raw(&M->vcenter[3*path[start]], axyz, arad, na)
                >= interior_threshold) break;
        if (start == np) return 0;
    }

    out[n++] = path[start];
    for (i = start + 1; i < np; i++) {
        double dx = M->center[3*path[i]]   - M->center[3*path[i-1]];
        double dy = M->center[3*path[i]+1] - M->center[3*path[i-1]+1];
        double dz = M->center[3*path[i]+2] - M->center[3*path[i-1]+2];
        if (sqrt(dx*dx+dy*dy+dz*dz) > 0.7) out[n++] = path[i];
    }
    while (n > 0) {
        int t = out[n-1];
        double dx = M->vcenter[3*t]-M->center[3*t];
        double dy = M->vcenter[3*t+1]-M->center[3*t+1];
        double dz = M->vcenter[3*t+2]-M->center[3*t+2];
        if (contains_point(M, t, &M->vcenter[3*t])
            || sqrt(dx*dx+dy*dy+dz*dz) < 3.0) break;
        n--;
    }
    return n < 5 ? 0 : n;
}

int mole_control_path(const mole_complex *M, const int *path, int np,
                      const double *axyz, const double *arad, int na,
                      double interior_threshold, int *out)
{
    return mole_control_path_ex(M, path, np, axyz, arad, na,
                                interior_threshold, out, 0);
}

int mole_profile(const mole_complex *M, const int *path, int np,
                 const double *axyz, const double *arad, int na,
                 mole_tunnel_profile *out)
{
    double *ts, *px, *py, *pz, *rs, *frs = NULL, *brs = NULL, d, len = 0.0;
    int i;
    const int NUM = 100;
    if (np < 2) return -1;

    ts = xa_malloc((size_t)np*sizeof(double));
    px = malloc((size_t)np*sizeof(double));
    py = malloc((size_t)np*sizeof(double));
    pz = malloc((size_t)np*sizeof(double));
    if (!ts||!px||!py||!pz) return -1;
    for (i = 0; i < np; i++) {
        ts[i] = (double)i / (np-1);
        px[i] = M->vcenter[3*path[i]];
        py[i] = M->vcenter[3*path[i]+1];
        pz[i] = M->vcenter[3*path[i]+2];
    }
    mole_spline_init(&out->sx, ts, px, np);
    mole_spline_init(&out->sy, ts, py, np);
    mole_spline_init(&out->sz, ts, pz, np);
    free(ts); free(px); free(py); free(pz);

    rs = xa_malloc((size_t)NUM*sizeof(double));
    ts = malloc((size_t)NUM*sizeof(double));
    frs = malloc((size_t)NUM*sizeof(double));
    brs = malloc((size_t)NUM*sizeof(double));
    if (!rs || !ts || !frs || !brs) return -1;
    d = 1.0/(NUM-1);
    for (i = 0; i < NUM; i++) {
        double t = d*i, u[3];
        ts[i] = t;
        u[0] = mole_spline_eval(&out->sx, t);
        u[1] = mole_spline_eval(&out->sy, t);
        u[2] = mole_spline_eval(&out->sz, t);
        rs[i] = radius_at(u, axyz, arad, na);
        frs[i] = free_radius_at(u, axyz, arad, na);
        brs[i] = bradius_at(u, axyz, arad, na);
        if (i > 0) {
            double v[3], t0 = d*(i-1);
            v[0] = mole_spline_eval(&out->sx, t0);
            v[1] = mole_spline_eval(&out->sy, t0);
            v[2] = mole_spline_eval(&out->sz, t0);
            len += sqrt((u[0]-v[0])*(u[0]-v[0]) + (u[1]-v[1])*(u[1]-v[1])
                      + (u[2]-v[2])*(u[2]-v[2]));
        }
    }
    mole_spline_init(&out->sr, ts, rs, NUM);
    mole_spline_init(&out->sfr, ts, frs, NUM);
    mole_spline_init(&out->sbr, ts, brs, NUM);
    out->length = len;
    free(rs); free(ts); free(frs); free(brs);
    return 0;
}

void mole_profile_free(mole_tunnel_profile *p)
{
    mole_spline_free(&p->sx); mole_spline_free(&p->sy);
    mole_spline_free(&p->sz); mole_spline_free(&p->sr);
    mole_spline_free(&p->sfr); mole_spline_free(&p->sbr);
}

/* ---- filters ----------------------------------------------------------- */

/* Complex.FilterBottleneck (TunnelComputation.cs:28-62).
 *
 * Skip the leading profile points whose DISTANCE is below BottleneckRadius -
 * MOLE reuses the radius parameter as a distance there, which is odd and is
 * theirs. If that skips everything the tunnel is KEPT (their `profile.Length ==
 * 0` returns true), not dropped.
 *
 * With tolerance 0 the tunnel survives only if nothing that remains is narrower
 * than BottleneckRadius. With a tolerance, a narrow stretch is allowed as long
 * as it is shorter than the tolerance; note MOLE seeds `start` at 0.0 rather
 * than at the first point's own distance, so the first stretch is measured from
 * the profile origin.
 */
int mole_filter_bottleneck(const mole_tunnel_profile *p, double bottleneck,
                           double tolerance, double density)
{
    /* GetProfile(d) yields (int)(length*d) + 1 samples at dt = 1/n, not n
       samples at 1/(n-1). Off by one row per profile against their CSVs. */
    int n = (int)(p->length * density), i;
    int started = 0, first = 1, have_start = 0;
    double dt, start = 0.0, dist = 0.0, px = 0, py = 0, pz = 0;
    /* n == 0 needs length < 1/density (0.5 A at the loosest gate, GetProfile(2)).
       MOLE's profile there is a single node, every field NaN (dt = 1/0 = +Inf,
       then dt*0 = NaN). Its own comparisons on that node are all "NaN < x",
       which IEEE makes false unconditionally - SkipWhile never skips it,
       Any() over it is false when tolerance is 0, and the tolerance branch's
       Skip(1) loop body never runs on a 1-element list. Every path returns
       true: MOLE keeps a fully-degenerate tunnel UNCONDITIONALLY. That is a
       derived verdict, matched directly rather than by attempting the NaN
       arithmetic itself - Tcl cannot: Inf*0 raises a domain error where IEEE
       yields NaN. The old clamp (n=1, then evaluate two real points) was
       WRONG: mole_filter_bottleneck_test.c shows it rejects 4 of 8 degenerate
       cases MOLE would keep. Not observed in the tested structures - shortest
       candidate over the six fixtures is 3.47 A against this 0.5 A gate, over
       247 candidates - not proven unreachable: the control path's minimum
       tetrahedron count is not a length floor, and adjacent sliver or
       near-degenerate tetrahedra can have arbitrarily close circumcentres.
       Watched by check 8b so a structure that ever approaches it says so. */
    if (n < 1) return 1;
    dt = 1.0/n;
    for (i = 0; i <= n; i++) {
        double t = dt*i, cx, cy, cz;
        double r = mole_spline_eval(&p->sr, t);
        /* TunnelProfile.Node.Distance is accumulated chord length
           (GetProfileNode: prev.Distance + prev.Center.DistanceTo(p)), not
           t * Length - the spline is not arc-length parameterised. SkipWhile(n
           => n.Distance < BottleneckRadius) keys on it. */
        cx = mole_spline_eval(&p->sx, t);
        cy = mole_spline_eval(&p->sy, t);
        cz = mole_spline_eval(&p->sz, t);
        if (i > 0) {
            double ax = cx-px, ay = cy-py, az = cz-pz;
            dist += sqrt(ax*ax + ay*ay + az*az);
        }
        px = cx; py = cy; pz = cz;
        if (!started) { if (dist < bottleneck) continue; started = 1; }
        if (tolerance == 0.0) {
            if (r < bottleneck) return 0;
            continue;
        }
        if (first) {
            /* profile.First(): seeds start at 0.0, not at its own distance. */
            first = 0;
            if (r < bottleneck) { have_start = 1; start = 0.0; }
            continue;
        }
        if (r < bottleneck) {
            if (have_start) {
                if (dist - start < tolerance) continue;
                return 0;
            }
            have_start = 1; start = dist;
        } else {
            have_start = 0;
        }
    }
    /* Everything skipped: MOLE keeps the tunnel. */
    if (!started && getenv("MOLE_BRANCH_AUDIT"))
        fprintf(stderr, "BRANCH empty-profile-kept len=%.3f\n", p->length);
    return 1;
}

/* Complex.FilterTunnels' similarity pass.
 *
 * Order by PATH length (vertex count), then for each pair take the one with the
 * shorter path and ask what fraction of its 6-per-Angstrom profile points lie
 * within 1 A of the other. Above MaxTunnelSimilarity the LONGER one is dropped.
 * Asymmetric on purpose: a long route that merely extends a short one is the
 * redundant member of the pair, not the other way round.
 */
static int sample_centreline(const mole_tunnel_profile *p, double density, double **out)
{
    int n = (int)(p->length * density), i;
    double dt;
    if (n < 1) n = 1;
    dt = 1.0 / n;
    *out = xa_malloc((size_t)(n + 1) * 3 * sizeof(double));
    if (!*out) return 0;
    for (i = 0; i <= n; i++) {
        double t = dt * i;
        (*out)[3*i+0] = mole_spline_eval(&p->sx, t);
        (*out)[3*i+1] = mole_spline_eval(&p->sy, t);
        (*out)[3*i+2] = mole_spline_eval(&p->sz, t);
    }
    return n + 1;
}

void mole_filter_similar(const mole_tunnel_profile *prof, const int *path_len,
                         char *dead, int n, double max_similarity)
{
    int i, j, *ord = malloc((size_t)(n ? n : 1) * sizeof(int));
    if (!ord) return;
    for (i = 0; i < n; i++) { dead[i] = 0; ord[i] = i; }

    /* Sort by PATH length ascending, then sweep i<j. A removed tunnel is
       skipped as pivot and as candidate, so the removal SEQUENCE decides which
       of a mutually-similar group survives. LINQ's OrderBy is stable, so ties
       hold their original order. */
    for (i = 1; i < n; i++) {
        int v = ord[i], k = i - 1;
        while (k >= 0 && path_len[ord[k]] > path_len[v]) { ord[k+1] = ord[k]; k--; }
        ord[k+1] = v;
    }

    for (i = 0; i < n; i++)
        for (j = i + 1; j < n; j++) {
            int a = ord[i], b = ord[j], na, nb, q, s, hit = 0;
            double *la, *lb;
            if (dead[a]) break;
            if (dead[b]) continue;
            /* RemoveLonger samples the shorter tunnel at GetProfile(6) and
               tests it against a tree built from the longer's GetProfile(2).
               The densities are asymmetric in MOLE. */
            na = sample_centreline(&prof[a], 6.0, &la);
            nb = sample_centreline(&prof[b], 2.0, &lb);
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
            if (getenv("MOLE_SIM_DEBUG"))
                fprintf(stderr, "  %.2f/p%d (pivot) vs %.2f/p%d : %d/%d = %.4f%s\n",
                        prof[a].length, path_len[a], prof[b].length, path_len[b],
                        hit, na, (double)hit/na,
                        (double)hit/na > max_similarity ? "  REMOVE" : "");
            if ((double)hit/na > max_similarity) dead[b] = 1;
            free(la); free(lb);
        }
    free(ord);
}
