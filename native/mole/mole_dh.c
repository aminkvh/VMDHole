/* MOLE 2's Delaunay triangulation, ported from WebChemistry.Framework.Core
 * (MIT): Geometry/Triangulation/DH/{DHTriangulation,HilbertOrdering,Tetrahedron,
 * DisconnectedFace,TriangulationVertex}.cs and DelaunayTriangulation.cs.
 *
 * Used instead of vor_delaunay because MOLE's per-cell VERTEX ORDER is part of
 * its output: the first control-path sample sits on a circumcentre equidistant
 * from its four atoms, so vertex order breaks that tie and decides a lining
 * layer boundary (Tunnel.cs:294 keys layers on an ordered identifier list).
 * The order is insertion history, so it is reproduced by reproducing the
 * insertion.
 *
 * Branch order and list discipline follow the C# because they are what the
 * history depends on. Pointers become pool indices; -1 is null.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_dh.h"
#include "../xalloc.h"

#define DH_NULL (-1)

/* ---------------------------------------------------------------- vertices */

/* TriangulationVertex. `index` is the ORIGINAL input position, assigned before
   the Hilbert reordering, and is what the face keys are built from. The single
   infinity vertex has Unit 0, LengthSquared 1 and index -1. */
typedef struct {
    double x, y, z, lensq, unit;
    int    index;
} dh_vert;

struct dh_face;

/* ------------------------------------------------------------ tetrahedra */

typedef struct {
    int    v[4], n[4];
    int    prev, next;          /* the C#'s Previous/Next, reused for the free list */
    int    opposite;
    double mu, mx, my, mz, mq;  /* MinorsUnit, MinorsX, MinorsY, MinorsZ, MinorsQ */
    double lastdet;
    unsigned char precomputed, infinite, marked, localflag;
    int    tag;
} dh_tet;

typedef struct {
    dh_vert *V;                 /* n real vertices, then INF at index n */
    int      n, inf;
    dh_tet  *T;
    int      ncap, nused, freetop;
    int     *simplices;
    int      nsimp, simpcap;
    int      faces_buf[1024];   /* FacesBuffer */
    struct dh_face *FA;         /* face pool; see below */
    int      fa_cap, fa_used, fa_free;
    int      failed;
} dh_ctx;

/* --------------------------------------------------------- disconnected faces */

/* DisconnectedFace: the three vertex INDICES of one face, sorted ascending, with
   the C#'s hash. Pooled the same way the tetrahedra are. */
struct dh_face {
    int tet, face;
    int v0, v1, v2;
    int hash;
    int prev;
};

/* ------------------------------------------------------------------ pools */

static int tet_alloc(dh_ctx *C)
{
    int t;
    if (C->freetop != DH_NULL) {
        t = C->freetop;
        C->freetop = C->T[t].prev;
        return t;
    }
    if (C->nused == C->ncap) {
        int cap = C->ncap ? C->ncap * 2 : 1024;
        dh_tet *p = realloc(C->T, (size_t)cap * sizeof(dh_tet));
        if (!p) { C->failed = 1; return DH_NULL; }
        C->T = p; C->ncap = cap;
    }
    return C->nused++;
}

/* TetrahedronFactory.Create / Tetrahedron.Recycle - one path, since a fresh
   slot and a recycled one are initialised identically here. */
static int tet_create(dh_ctx *C, int v0, int v1, int v2, int v3, int opposite)
{
    int t = tet_alloc(C);
    dh_tet *p;
    if (t == DH_NULL) return DH_NULL;
    p = &C->T[t];
    p->v[0] = v0; p->v[1] = v1; p->v[2] = v2; p->v[3] = v3;
    p->infinite = (unsigned char)(C->V[v0].unit == 0 || C->V[v1].unit == 0 ||
                                  C->V[v2].unit == 0 || C->V[v3].unit == 0);
    p->n[0] = p->n[1] = p->n[2] = p->n[3] = DH_NULL;
    p->opposite = opposite;
    p->precomputed = p->marked = p->localflag = 0;
    p->prev = p->next = DH_NULL;
    p->lastdet = 0.0;
    p->tag = -1;
    return t;
}

static void tet_dispose(dh_ctx *C, int t)
{
    C->T[t].prev = C->freetop;
    C->freetop = t;
}

/* Tetrahedron.VertexUpdated */
static void tet_vertex_updated(dh_ctx *C, int t)
{
    dh_tet *p = &C->T[t];
    p->infinite = (unsigned char)(C->V[p->v[0]].unit == 0 || C->V[p->v[1]].unit == 0 ||
                                  C->V[p->v[2]].unit == 0 || C->V[p->v[3]].unit == 0);
    p->precomputed = 0;
}

/* Tetrahedron.UpdateLink. Note the C#'s final `else` is unconditional: a link
   that matches none of N0..N2 is assumed to be N3. */
static void tet_update_link(dh_ctx *C, int t, int oldl, int newl)
{
    dh_tet *p = &C->T[t];
    if      (p->n[0] == oldl) p->n[0] = newl;
    else if (p->n[1] == oldl) p->n[1] = newl;
    else if (p->n[2] == oldl) p->n[2] = newl;
    else                      p->n[3] = newl;
}

/* Tetrahedron.Precompute, in MOLE's operation order: six shared 2x2 minors per
   vertex pair, five 4x4 cofactors, two sign normalisations. The order fixes the
   last bits, which decide the vertex-order ties this file reproduces. */
static void tet_precompute(dh_ctx *C, int t)
{
    dh_tet *p = &C->T[t];
    const dh_vert *V0 = &C->V[p->v[0]], *V1 = &C->V[p->v[1]];
    const dh_vert *V2 = &C->V[p->v[2]], *V3 = &C->V[p->v[3]];
    double m00, m01, m02, m03, m04, m05, m06, m07, m08, m09;
    double m10, m11, m12, m13, m14, m15, m16, m17, m18, m19;

    p->precomputed = 1;

    m00 = V0->unit * V1->x     - V0->x     * V1->unit;
    m01 = V0->unit * V1->y     - V0->y     * V1->unit;
    m02 = V0->unit * V1->z     - V0->z     * V1->unit;
    m03 = V0->unit * V1->lensq - V0->lensq * V1->unit;
    m04 = V0->x    * V1->y     - V0->y     * V1->x;
    m05 = V0->x    * V1->z     - V0->z     * V1->x;
    m06 = V0->x    * V1->lensq - V0->lensq * V1->x;
    m07 = V0->y    * V1->z     - V0->z     * V1->y;
    m08 = V0->y    * V1->lensq - V0->lensq * V1->y;
    m09 = V0->z    * V1->lensq - V0->lensq * V1->z;

    m10 = V2->unit * V3->x     - V2->x     * V3->unit;
    m11 = V2->unit * V3->y     - V2->y     * V3->unit;
    m12 = V2->unit * V3->z     - V2->z     * V3->unit;
    m13 = V2->unit * V3->lensq - V2->lensq * V3->unit;
    m14 = V2->x    * V3->y     - V2->y     * V3->x;
    m15 = V2->x    * V3->z     - V2->z     * V3->x;
    m16 = V2->x    * V3->lensq - V2->lensq * V3->x;
    m17 = V2->y    * V3->z     - V2->z     * V3->y;
    m18 = V2->y    * V3->lensq - V2->lensq * V3->y;
    m19 = V2->z    * V3->lensq - V2->lensq * V3->z;

    p->mu = m04*m19 - m05*m18 + m06*m17 + m07*m16 - m08*m15 + m09*m14;
    p->mx = m01*m19 - m02*m18 + m03*m17 + m07*m13 - m08*m12 + m09*m11;
    p->my = m00*m19 - m02*m16 + m03*m15 + m05*m13 - m06*m12 + m09*m10;
    p->mz = m00*m18 - m01*m16 + m03*m14 + m04*m13 - m06*m11 + m08*m10;
    p->mq = m00*m17 - m01*m15 + m02*m14 + m04*m12 - m05*m11 + m07*m10;

    if (p->mq < 0) {
        p->mu = -p->mu; p->mx = -p->mx; p->my = -p->my;
        p->mz = -p->mz; p->mq = -p->mq;
    } else if (p->infinite) {
        const dh_vert *O = &C->V[p->opposite];
        if (p->mu - p->mx*O->x + p->my*O->y - p->mz*O->z < 0) {
            p->mu = -p->mu; p->mx = -p->mx; p->my = -p->my; p->mz = -p->mz;
        }
    }
}

static double tet_sphere_det(dh_ctx *C, int t, int vtx)
{
    dh_tet *p = &C->T[t];
    const dh_vert *v = &C->V[vtx];
    if (!p->precomputed) tet_precompute(C, t);
    return p->mu - p->mx*v->x + p->my*v->y - p->mz*v->z + p->mq*v->lensq;
}

/* ---------------------------------------------------------- Hilbert ordering */

/* Hamilton, "Compact Hilbert Indices" (2006) p.19, 3D specialisation. MOLE
   ships tables for orders 1-4 and falls back to this above; the tables are
   exactly what this produces, so it is used for every order. */
int dh_hilbert_encode(int order, int x, int y, int z)
{
    static const unsigned bitMask[5] = { 0x1, 0x2, 0x4, 0x8, 0x10 };
    static const unsigned E[8] = { 0, 0, 0, 3, 3, 6, 6, 5 };
    static const unsigned D[8] = { 0, 1, 1, 2, 2, 1, 1, 0 };
    unsigned h = 0, e = 0, l, w;
    int d = 1, i;
    for (i = order - 1; i >= 0; i--) {
        l  = ((unsigned)x & bitMask[i]) >> i;
        l |= (((unsigned)y & bitMask[i]) >> i) << 1;
        l |= (((unsigned)z & bitMask[i]) >> i) << 2;
        l ^= e;
        l = (unsigned)((l >> (d+1)) | (l << (3-(d+1)))) & 0x7;   /* RotRight3b */
        w = l ^ (l >> 1) ^ (l >> 2);                             /* GrayCodeInverse3b */
        w &= 0x7;
        e ^= (unsigned)((E[w] << (d+1)) | (E[w] >> (3-(d+1)))) & 0x7;  /* RotLeft3b */
        d = (int)((d + D[w] + 1) % 3);
        h = (h << 3) | w;
    }
    return (int)h;
}

#define DH_VOXEL_THRESHOLD 400

static int dh_curve_order(int n)
{
    if (n < DH_VOXEL_THRESHOLD)          return 0;
    else if (n < DH_VOXEL_THRESHOLD*8)   return 1;
    else if (n < DH_VOXEL_THRESHOLD*64)  return 2;
    else if (n < DH_VOXEL_THRESHOLD*512) return 3;
    else                                 return 4;
}

/* HilbertOrdering.OrderPointsSegment. `src` holds vertex slots; the ordered
   result is written into target[start..start+count-1]. */
static int order_segment(const double *xyz, const int *src, int *target, int start,
                         double minx, double miny, double minz,
                         double maxx, double maxy, double maxz, int count)
{
    int order = dh_curve_order(count), cells, nvox, i, j, rc = 0;
    int *hidx = NULL, *dens = NULL, *cum = NULL, *lut = NULL, *newsrc = NULL;
    double stepx, stepy, stepz;

    if (order == 0) {
        for (i = 0; i < count; i++) target[start + i] = src[i];
        return 0;
    }
    cells = 1 << order;
    nvox  = 1 << (3 * order);
    hidx = xa_malloc((size_t)count * sizeof(int));
    dens = calloc((size_t)nvox, sizeof(int));
    cum  = calloc((size_t)nvox, sizeof(int));
    lut  = malloc((size_t)nvox * 3 * sizeof(int));
    if (!hidx || !dens || !cum || !lut) { rc = -1; goto done; }

    for (j = 0; j < count; j++) {
        const double *p = &xyz[3 * src[j]];
        /* ConvertCoordinates: C#'s (int) cast truncates toward zero, and max was
           inflated by 0.001*(max-min) so the quotient stays below 1. */
        int cx = (int)(((p[0] - minx) / (maxx - minx)) * cells);
        int cy = (int)(((p[1] - miny) / (maxy - miny)) * cells);
        int cz = (int)(((p[2] - minz) / (maxz - minz)) * cells);
        int h  = dh_hilbert_encode(order, cx, cy, cz);
        dens[h]++;
        lut[3*h] = cx; lut[3*h+1] = cy; lut[3*h+2] = cz;
        hidx[j] = h;
    }
    for (i = 1; i < nvox; i++) cum[i] = cum[i-1] + dens[i-1];
    for (j = 0; j < count; j++) {
        target[start + cum[hidx[j]]] = src[j];
        cum[hidx[j]]++;
    }

    stepx = (maxx - minx) / cells;
    stepy = (maxy - miny) / cells;
    stepz = (maxz - minz) / cells;
    /* From 1, not 0: MOLE never re-orders voxel 0 however crowded it is. */
    for (i = 1; i < nvox; i++) {
        if (dens[i] >= DH_VOXEL_THRESHOLD) {
            int base = start + cum[i] - dens[i];
            double lx, ly, lz, hx, hy, hz;
            newsrc = malloc((size_t)dens[i] * sizeof(int));
            if (!newsrc) { rc = -1; goto done; }
            memcpy(newsrc, &target[base], (size_t)dens[i] * sizeof(int));
            lx = lut[3*i]  *stepx + minx; hx = (lut[3*i]  +1)*stepx + minx;
            ly = lut[3*i+1]*stepy + miny; hy = (lut[3*i+1]+1)*stepy + miny;
            lz = lut[3*i+2]*stepz + minz; hz = (lut[3*i+2]+1)*stepz + minz;
            rc = order_segment(xyz, newsrc, target, base, lx, ly, lz, hx, hy, hz, dens[i]);
            free(newsrc); newsrc = NULL;
            if (rc) goto done;
        }
    }
done:
    free(hidx); free(dens); free(cum); free(lut); free(newsrc);
    return rc;
}

/* HilbertOrdering.OrderPoints */
static int order_points(const double *xyz, int n, int *out)
{
    double minx = 1e308, miny = 1e308, minz = 1e308;
    double maxx = -1e308, maxy = -1e308, maxz = -1e308;
    int *src = malloc((size_t)n * sizeof(int));
    int i, rc;
    if (!src) return -1;
    for (i = 0; i < n; i++) {
        const double *p = &xyz[3*i];
        if (p[0] > maxx) maxx = p[0];
        if (p[1] > maxy) maxy = p[1];
        if (p[2] > maxz) maxz = p[2];
        if (p[0] < minx) minx = p[0];
        if (p[1] < miny) miny = p[1];
        if (p[2] < minz) minz = p[2];
        src[i] = i;
    }
    maxx += 0.001 * (maxx - minx);
    maxy += 0.001 * (maxy - miny);
    maxz += 0.001 * (maxz - minz);
    rc = order_segment(xyz, src, out, 0, minx, miny, minz, maxx, maxy, maxz, n);
    free(src);
    return rc;
}

/* ---------------------------------------------------------- face bookkeeping */

static int face_create(dh_ctx *C, int tet, int face)
{
    struct dh_face *FA;
    int f, a, b, c, t, hash;
    dh_tet *p = &C->T[tet];
    if (C->fa_free != DH_NULL) { f = C->fa_free; C->fa_free = C->FA[f].prev; }
    else {
        if (C->fa_used == C->fa_cap) {
            int cap = C->fa_cap ? C->fa_cap * 2 : 1024;
            struct dh_face *q = realloc(C->FA, (size_t)cap * sizeof(struct dh_face));
            if (!q) { C->failed = 1; return DH_NULL; }
            C->FA = q; C->fa_cap = cap;
        }
        f = C->fa_used++;
    }
    FA = C->FA;
    switch (face) {
        case 0: a = C->V[p->v[1]].index; b = C->V[p->v[2]].index; c = C->V[p->v[3]].index; break;
        case 1: a = C->V[p->v[0]].index; b = C->V[p->v[2]].index; c = C->V[p->v[3]].index; break;
        case 2: a = C->V[p->v[0]].index; b = C->V[p->v[1]].index; c = C->V[p->v[3]].index; break;
        default:a = C->V[p->v[0]].index; b = C->V[p->v[1]].index; c = C->V[p->v[2]].index; break;
    }
    /* MOLE's own three-swap sort, kept as written. */
    if (b < a) { t = b; b = a; a = t; }
    if (c < b) { t = c; c = b; b = t; }
    if (b < a) { t = b; b = a; a = t; }
    hash = 23;
    hash = 31 * hash + a;
    hash = 31 * hash + b;
    hash = 31 * hash + c;
    FA[f].tet = tet; FA[f].face = face;
    FA[f].v0 = a; FA[f].v1 = b; FA[f].v2 = c;
    FA[f].hash = hash;
    FA[f].prev = DH_NULL;
    return f;
}

static void face_dispose(dh_ctx *C, int f) { C->FA[f].prev = C->fa_free; C->fa_free = f; }

static int add_new_face(dh_ctx *C, int list, int tet, int face)
{
    int f = face_create(C, tet, face);
    if (f == DH_NULL) return list;
    C->FA[f].prev = list;
    return f;
}

/* ConnectFaces: face1.Face*4 + face2.Face selects which slot on each side. */
static void connect_faces(dh_ctx *C, int f1, int f2)
{
    struct dh_face *FA = C->FA;
    C->T[FA[f2].tet].n[FA[f2].face] = FA[f1].tet;
    C->T[FA[f1].tet].n[FA[f1].face] = FA[f2].tet;
}

static int can_connect(dh_ctx *C, int f1, int f2)
{
    struct dh_face *FA = C->FA;
    return FA[f1].v0 == FA[f2].v0 && FA[f1].v1 == FA[f2].v1 && FA[f1].v2 == FA[f2].v2;
}

/* DHTriangulation.CreateLinks */
static void create_links(dh_ctx *C, int faces)
{
    struct dh_face *FA = C->FA;
    while (faces != DH_NULL) {
        int cur = faces, hash;
        faces = FA[faces].prev;
        hash = FA[cur].hash & 0x3FF;

        if (C->faces_buf[hash] == DH_NULL) {
            C->faces_buf[hash] = cur;
            FA[cur].prev = DH_NULL;
        } else if (FA[C->faces_buf[hash]].hash == FA[cur].hash &&
                   can_connect(C, C->faces_buf[hash], cur)) {
            int tmp = C->faces_buf[hash];
            C->faces_buf[hash] = FA[tmp].prev;
            connect_faces(C, tmp, cur);
            face_dispose(C, tmp);
            face_dispose(C, cur);
        } else {
            int last = C->faces_buf[hash], f, hit = 0;
            for (f = FA[last].prev; f != DH_NULL; f = FA[f].prev) {
                if (FA[cur].hash == FA[f].hash && can_connect(C, f, cur)) {
                    connect_faces(C, f, cur);
                    FA[last].prev = FA[f].prev;
                    face_dispose(C, f);
                    face_dispose(C, cur);
                    hit = 1;
                    break;
                }
                last = f;
            }
            if (!hit) {
                FA[cur].prev = C->faces_buf[hash];
                C->faces_buf[hash] = cur;
            }
        }
    }
}

/* ------------------------------------------------------------------ walking */

/* FindNextStepDegenerated's shared 2x2 minors, and the per-face determinant. */
static double plane_det(dh_ctx *C, int t, int face, const dh_vert *vx, double *qout)
{
    dh_tet *p = &C->T[t];
    const dh_vert *V0 = &C->V[p->v[0]], *V1 = &C->V[p->v[1]];
    const dh_vert *V2 = &C->V[p->v[2]], *V3 = &C->V[p->v[3]];
    double m00 = V1->x - V0->x, m01 = V1->y - V0->y, m02 = V1->z - V0->z;
    double m03 = V0->x*V1->y - V1->x*V0->y;
    double m04 = V0->x*V1->z - V1->x*V0->z;
    double m05 = V0->y*V1->z - V1->y*V0->z;
    double m10 = V3->x - V2->x, m11 = V3->y - V2->y, m12 = V3->z - V2->z;
    double m13 = V2->x*V3->y - V3->x*V2->y;
    double m14 = V2->x*V3->z - V3->x*V2->z;
    double m15 = V2->y*V3->z - V3->y*V2->z;
    if (qout) *qout = m00*m15 - m01*m14 + m02*m13 + m03*m12 - m04*m11 + m05*m10;
    switch (face) {
        case 0: return (V1->x - vx->x)*m15 - (V1->y - vx->y)*m14 + (V1->z - vx->z)*m13
                     + (vx->x*V1->y - V1->x*vx->y)*m12 - (vx->x*V1->z - V1->x*vx->z)*m11
                     + (vx->y*V1->z - V1->y*vx->z)*m10;
        case 1: return (vx->x - V0->x)*m15 - (vx->y - V0->y)*m14 + (vx->z - V0->z)*m13
                     + (V0->x*vx->y - vx->x*V0->y)*m12 - (V0->x*vx->z - vx->x*V0->z)*m11
                     + (V0->y*vx->z - vx->y*V0->z)*m10;
        case 2: return m00*(vx->y*V3->z - V3->y*vx->z) - m01*(vx->x*V3->z - V3->x*vx->z)
                     + m02*(vx->x*V3->y - V3->x*vx->y) + m03*(V3->z - vx->z)
                     - m04*(V3->y - vx->y) + m05*(V3->x - vx->x);
        default:return m00*(V2->y*vx->z - vx->y*V2->z) - m01*(V2->x*vx->z - vx->x*V2->z)
                     + m02*(V2->x*vx->y - vx->x*V2->y) + m03*(vx->z - V2->z)
                     - m04*(vx->y - V2->y) + m05*(vx->x - V2->x);
    }
}

static int find_next_step_degenerated(dh_ctx *C, int cur, int prev, double *det, int vtx)
{
    const dh_vert *vx = &C->V[vtx];
    double q, d;
    int k;
    for (k = 0; k < 4; k++) {
        int nb = C->T[cur].n[k];
        if (nb == prev) continue;
        d = plane_det(C, cur, k, vx, &q);
        if (d * q < 0) {
            *det = tet_sphere_det(C, nb, vtx);
            return nb;
        }
    }
    /* MOLE throws InvalidOperationException here. */
    C->failed = 1;
    return cur;
}

/* PlaneCompareDegenerated */
static int plane_compare_degenerated(dh_ctx *C, int t, int face, int vtx)
{
    double q, d = plane_det(C, t, face, &C->V[vtx], &q);
    return d * q < 0;
}

static int find_next_step(dh_ctx *C, int cur, int prev, double *det, int vtx)
{
    int k;
    for (k = 0; k < 4; k++) {
        int nb = C->T[cur].n[k];
        double newdet;
        if (nb == prev) continue;
        newdet = tet_sphere_det(C, nb, vtx);
        if (C->T[nb].mq * (*det) - C->T[cur].mq * newdet > 0) {
            *det = newdet;
            return nb;
        }
    }
    return find_next_step_degenerated(C, cur, prev, det, vtx);
}

/* DHTriangulation.SearchNext */
static void search_next(dh_ctx *C, int nb, int *stack)
{
    if (!C->T[nb].marked) {
        C->T[nb].marked = 1;
        C->T[nb].prev = *stack;
        *stack = nb;
    }
}

/* DHTriangulation.Localize */
static int localize(dh_ctx *C, int vtx, int cur)
{
    int previous = DH_NULL, pathlen = 0, tosearch, neighbours, result = DH_NULL, t;
    double lastdet = tet_sphere_det(C, cur, vtx);

    while (lastdet > 0) {
        int tmp = cur;
        cur = (pathlen < 500) ? find_next_step(C, cur, previous, &lastdet, vtx)
                              : find_next_step_degenerated(C, cur, previous, &lastdet, vtx);
        previous = tmp;
        pathlen++;
        if (C->failed) return DH_NULL;
    }

    neighbours = previous;
    if (neighbours != DH_NULL) {
        C->T[neighbours].prev = DH_NULL;
        C->T[neighbours].marked = 1;
        C->T[neighbours].lastdet = tet_sphere_det(C, neighbours, vtx);
    }
    tosearch = cur;
    C->T[tosearch].prev = DH_NULL;
    C->T[tosearch].marked = 1;

    while (tosearch != DH_NULL) {
        cur = tosearch;
        tosearch = C->T[cur].prev;
        C->T[cur].lastdet = tet_sphere_det(C, cur, vtx);
        if (C->T[cur].lastdet <= 0) {
            if (result != DH_NULL) C->T[result].next = cur;
            C->T[cur].prev = result;
            result = cur;
            search_next(C, C->T[cur].n[0], &tosearch);
            search_next(C, C->T[cur].n[1], &tosearch);
            search_next(C, C->T[cur].n[2], &tosearch);
            search_next(C, C->T[cur].n[3], &tosearch);
        } else {
            C->T[cur].prev = neighbours;
            neighbours = cur;
        }
    }

    for (t = neighbours; t != DH_NULL; t = C->T[t].prev) C->T[t].marked = 0;

    /* "Check if the cavity is properly shaped". MOLE throws here; record the
       failure and unwind rather than continue on a broken triangulation. */
    for (t = neighbours; t != DH_NULL; t = C->T[t].prev) {
        int k;
        for (k = 0; k < 4; k++) {
            int nb = C->T[t].n[k];
            double lhs;
            if (nb == DH_NULL || !C->T[nb].marked || C->T[nb].infinite) continue;
            lhs = C->T[t].mq * C->T[nb].lastdet - C->T[nb].mq * C->T[t].lastdet;
            if (lhs > 1e-16) { C->failed = 2; return DH_NULL; }
            if (lhs > -1e-16 && !plane_compare_degenerated(C, t, k, vtx)) {
                C->failed = 2; return DH_NULL;
            }
        }
    }
    return result;
}

/* DHTriangulation.RemoveFromList */
static void remove_from_list(dh_ctx *C, int t, int *list)
{
    if (C->T[t].next != DH_NULL) {
        C->T[C->T[t].next].prev = C->T[t].prev;
        if (C->T[t].prev != DH_NULL) C->T[C->T[t].prev].next = C->T[t].next;
    } else {
        *list = C->T[t].prev;
        if (*list != DH_NULL) C->T[*list].next = DH_NULL;
    }
}

/* DHTriangulation.TryUpdate */
static int try_update(dh_ctx *C, int t, int opposite, int newv)
{
    dh_tet *p = &C->T[t];
    int slot;
    if      (p->v[0] == opposite) slot = 0;
    else if (p->v[1] == opposite) slot = 1;
    else if (p->v[2] == opposite) slot = 2;
    else                          slot = 3;
    if (p->n[slot] == DH_NULL || C->T[p->n[slot]].marked) return 0;
    p->v[slot] = newv;
    tet_vertex_updated(C, t);
    return 1;
}

/* The per-slot body of Insert's inner loop, which the C# writes out four times.
   The three vertices handed to Create are the current cell's OTHER three in
   ascending slot order, which is why the argument list differs per slot. */
static int insert_slot(dh_ctx *C, int cur, int slot, int oldv, int vtx,
                       int *toprocess, int *localstack, int *result, int *innerfaces)
{
    static const int OTHER[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
    dh_tet *p = &C->T[cur];
    int nb = p->n[slot];

    if (nb == DH_NULL) { *innerfaces = add_new_face(C, *innerfaces, cur, slot); return 0; }
    if (C->T[nb].localflag) return 0;

    if (!C->T[nb].marked) {
        const int *o = OTHER[slot];
        int fresh = tet_create(C, p->v[o[0]], p->v[o[1]], p->v[o[2]], oldv, p->v[slot]);
        if (fresh == DH_NULL) return -1;
        p = &C->T[cur];                       /* tet_create may have realloc'd */
        p->n[slot] = fresh;
        C->T[fresh].n[3] = cur;
        tet_update_link(C, nb, cur, fresh);
        if (*result == DH_NULL && !C->T[fresh].infinite) *result = fresh;
        if (C->T[fresh].v[0] == vtx) {
            C->T[fresh].n[0] = nb;
            *innerfaces = add_new_face(C, *innerfaces, fresh, 1);
            *innerfaces = add_new_face(C, *innerfaces, fresh, 2);
        } else if (C->T[fresh].v[1] == vtx) {
            C->T[fresh].n[1] = nb;
            *innerfaces = add_new_face(C, *innerfaces, fresh, 0);
            *innerfaces = add_new_face(C, *innerfaces, fresh, 2);
        } else if (C->T[fresh].v[2] == vtx) {
            C->T[fresh].n[2] = nb;
            *innerfaces = add_new_face(C, *innerfaces, fresh, 0);
            *innerfaces = add_new_face(C, *innerfaces, fresh, 1);
        }
    } else if (try_update(C, nb, oldv, vtx)) {
        remove_from_list(C, nb, toprocess);
        C->T[nb].prev = *localstack;
        *localstack = nb;
        C->T[nb].localflag = 1;
    } else {
        tet_update_link(C, nb, cur, DH_NULL);
        *innerfaces = add_new_face(C, *innerfaces, cur, slot);
        C->T[cur].n[slot] = DH_NULL;
    }
    return 0;
}

/* DHTriangulation.Insert */
static int insert(dh_ctx *C, int toprocess, int vtx)
{
    int innerfaces = DH_NULL, result = DH_NULL, tounmark = DH_NULL, t;

    while (toprocess != DH_NULL) {
        int cur = toprocess, localstack, oldv, slot;
        toprocess = C->T[cur].prev;
        if (toprocess != DH_NULL) C->T[toprocess].next = DH_NULL;

        slot = -1;
        for (t = 0; t < 4; t++) {
            int nb = C->T[cur].n[t];
            if (nb != DH_NULL && !C->T[nb].marked) { slot = t; break; }
        }
        if (slot < 0) {
            /* Every neighbour is inside the cavity - unlink and dispose. */
            for (t = 0; t < 4; t++)
                if (C->T[cur].n[t] != DH_NULL)
                    tet_update_link(C, C->T[cur].n[t], cur, DH_NULL);
            tet_dispose(C, cur);
            continue;
        }
        oldv = C->T[cur].v[slot];
        C->T[cur].v[slot] = vtx;
        tet_vertex_updated(C, cur);
        C->T[cur].prev = DH_NULL;
        localstack = cur;
        C->T[cur].localflag = 1;

        while (localstack != DH_NULL) {
            cur = localstack;
            localstack = C->T[cur].prev;
            if (result == DH_NULL && !C->T[cur].infinite) result = cur;
            C->T[cur].prev = tounmark;
            tounmark = cur;

            for (t = 0; t < 4; t++) {
                if (C->T[cur].v[t] == vtx) continue;
                if (insert_slot(C, cur, t, oldv, vtx, &toprocess, &localstack,
                                &result, &innerfaces) != 0) return DH_NULL;
            }
        }
    }

    for (t = tounmark; t != DH_NULL; t = C->T[t].prev) {
        C->T[t].marked = 0;
        C->T[t].localflag = 0;
    }
    create_links(C, innerfaces);
    return result;
}

/* DHTriangulation.Init */
static int dh_init(dh_ctx *C, int *ordered, int n)
{
    int first, c1, c2, c3, c4, i;

    first = tet_create(C, ordered[0], ordered[1], ordered[2], ordered[3], C->inf);
    if (first == DH_NULL) return -1;
    tet_sphere_det(C, first, C->inf);

    /* Four coplanar starting points. MOLE swaps in a random later vertex
       (`new Random()`), so its output is not reproducible when this fires. */
    i = 4;
    while (C->T[first].mq < 0.0000001) {
        int tmp;
        if (i >= n) return -1;
        fprintf(stderr, "mole_dh: first four points are coplanar; MOLE picks the "
                        "replacement at random here, so this run cannot be "
                        "compared with it\n");
        tmp = ordered[0]; ordered[0] = ordered[i]; ordered[i] = tmp;
        i++;
        tet_dispose(C, first);
        first = tet_create(C, ordered[0], ordered[1], ordered[2], ordered[3], C->inf);
        if (first == DH_NULL) return -1;
        tet_sphere_det(C, first, C->inf);
    }

    {
        int v0 = C->T[first].v[0], v1 = C->T[first].v[1];
        int v2 = C->T[first].v[2], v3 = C->T[first].v[3];
        c1 = tet_create(C, v0, v1, v2, C->inf, v3);
        c2 = tet_create(C, v1, v2, v3, C->inf, v0);
        c3 = tet_create(C, v2, v3, v0, C->inf, v1);
        c4 = tet_create(C, v3, v0, v1, C->inf, v2);
        if (c1 == DH_NULL || c2 == DH_NULL || c3 == DH_NULL || c4 == DH_NULL)
            return -1;
    }
    C->T[c1].n[3] = first; C->T[c2].n[3] = first;
    C->T[c3].n[3] = first; C->T[c4].n[3] = first;
    C->T[first].n[3] = c1; C->T[first].n[0] = c2;
    C->T[first].n[1] = c3; C->T[first].n[2] = c4;
    C->T[c1].n[2] = c4; C->T[c1].n[1] = c3; C->T[c1].n[0] = c2;
    C->T[c2].n[2] = c1; C->T[c2].n[1] = c4; C->T[c2].n[0] = c3;
    C->T[c3].n[2] = c2; C->T[c3].n[1] = c1; C->T[c3].n[0] = c4;
    C->T[c4].n[2] = c3; C->T[c4].n[1] = c2; C->T[c4].n[0] = c1;
    return first;
}

/* DHTriangulation.Finish: depth-first from `start`, which is what fixes the cell
   ORDER of the result and therefore our tetrahedron indices. */
static int dh_finish(dh_ctx *C, int start)
{
    int stack = start, k;
    C->T[start].prev = DH_NULL;
    C->T[start].marked = 1;
    while (stack != DH_NULL) {
        int t = stack;
        stack = C->T[t].prev;
        if (C->nsimp == C->simpcap) {
            int cap = C->simpcap ? C->simpcap * 2 : 1024;
            int *p = realloc(C->simplices, (size_t)cap * sizeof(int));
            if (!p) return -1;
            C->simplices = p; C->simpcap = cap;
        }
        C->simplices[C->nsimp++] = t;
        for (k = 0; k < 4; k++) {
            int nb = C->T[t].n[k];
            if (nb == DH_NULL) continue;
            if (!C->T[nb].marked) {
                C->T[nb].prev = stack;
                stack = nb;
                C->T[nb].marked = 1;
            }
        }
    }
    return 0;
}

/* ------------------------------------------------------------------- driver */

int dh_build(dh_mesh *out, const double *xyz, int n)
{
    dh_ctx C;
    int *ordered = NULL, i, start, rc = -1, nfin = 0;

    memset(out, 0, sizeof(*out));
    if (n < 4) return -1;
    memset(&C, 0, sizeof C);
    for (i = 0; i < 1024; i++) C.faces_buf[i] = DH_NULL;
    C.freetop = DH_NULL;
    C.fa_free = DH_NULL;

    C.n = n; C.inf = n;
    C.V = xa_malloc((size_t)(n + 1) * sizeof(dh_vert));
    ordered = malloc((size_t)n * sizeof(int));
    if (!C.V || !ordered) goto done;
    for (i = 0; i < n; i++) {
        C.V[i].x = xyz[3*i]; C.V[i].y = xyz[3*i+1]; C.V[i].z = xyz[3*i+2];
        C.V[i].unit = 1.0;
        C.V[i].lensq = C.V[i].x*C.V[i].x + C.V[i].y*C.V[i].y + C.V[i].z*C.V[i].z;
        C.V[i].index = i;
    }
    /* The infinity vertex: Unit 0, LengthSquared 1, coordinates 0, index -1. */
    C.V[n].x = C.V[n].y = C.V[n].z = 0.0;
    C.V[n].unit = 0.0; C.V[n].lensq = 1.0; C.V[n].index = -1;

    if (order_points(xyz, n, ordered) != 0) goto done;

    start = dh_init(&C, ordered, n);
    if (start == DH_NULL) goto done;
    for (i = 4; i < n; i++) {
        int region = localize(&C, ordered[i], start);
        if (C.failed) goto done;
        start = insert(&C, region, ordered[i]);
        if (C.failed || start == DH_NULL) goto done;
    }
    if (dh_finish(&C, start) != 0) goto done;

    /* DelaunayTriangulation3D.Create: finite cells only, tagged in list order,
       then adjacency mapped through the tags. */
    for (i = 0; i < C.nsimp; i++)
        if (!C.T[C.simplices[i]].infinite) C.T[C.simplices[i]].tag = nfin++;
    out->nt = nfin;
    out->tv = xa_malloc((size_t)(nfin ? nfin : 1) * 4 * sizeof(int));
    out->tn = xa_malloc((size_t)(nfin ? nfin : 1) * 4 * sizeof(int));
    if (!out->tv || !out->tn) goto done;
    for (i = 0; i < C.nsimp; i++) {
        int t = C.simplices[i], k, c = C.T[t].tag;
        if (C.T[t].infinite) continue;
        for (k = 0; k < 4; k++) {
            out->tv[4*c+k] = C.T[t].v[k];
            out->tn[4*c+k] = (C.T[t].n[k] != DH_NULL && !C.T[C.T[t].n[k]].infinite)
                             ? C.T[C.T[t].n[k]].tag : -1;
        }
    }
    rc = 0;
done:
    if (rc != 0 && C.failed == 2)
        fprintf(stderr, "mole_dh: degenerate cavity - MOLE throws here too\n");
    free(C.V); free(C.T); free(C.simplices); free(ordered); free(C.FA);
    if (rc != 0) { free(out->tv); free(out->tn); memset(out, 0, sizeof(*out)); }
    return rc;
}

void dh_mesh_free(dh_mesh *m)
{
    free(m->tv); free(m->tn);
    memset(m, 0, sizeof(*m));
}
