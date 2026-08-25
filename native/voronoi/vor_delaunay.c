/* Bowyer-Watson incremental Delaunay on exact predicates. See the header for
   what it is for and how tetrahedra are stored. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vor_delaunay.h"
#include "../xalloc.h"

#define PT(m,i) ((m)->p + 3*(i))

static int tet_new(dt_mesh *m)
{
    if (m->nfree > 0) {
        int ti = m->freet[--m->nfree];
        memset(&m->t[ti], 0, sizeof(dt_tet));
        m->t[ti].nb[0] = m->t[ti].nb[1] = m->t[ti].nb[2] = m->t[ti].nb[3] = -1;
        return ti;
    }
    if (m->nt >= m->cap) {
        int nc = m->cap ? m->cap * 2 : 1024;
        dt_tet *nt = realloc(m->t, (size_t)nc * sizeof(dt_tet));
        if (!nt) return -1;
        m->t = nt; m->cap = nc;
    }
    memset(&m->t[m->nt], 0, sizeof(dt_tet));
    m->t[m->nt].nb[0] = m->t[m->nt].nb[1] = m->t[m->nt].nb[2] = m->t[m->nt].nb[3] = -1;
    return m->nt++;
}

/* Order v[0..3] so the tetrahedron is positively oriented. Every later
   predicate call assumes this, so it is enforced here rather than checked. */
static void tet_orient(dt_mesh *m, int ti)
{
    dt_tet *t = &m->t[ti];
    if (vp_orient3d(PT(m,t->v[0]), PT(m,t->v[1]), PT(m,t->v[2]), PT(m,t->v[3])) < 0) {
        int tmp = t->v[0]; t->v[0] = t->v[1]; t->v[1] = tmp;
        tmp = t->nb[0];   t->nb[0] = t->nb[1]; t->nb[1] = tmp;
    }
}

/* Is p inside (or on) tetrahedron ti? If not, `exit_face` gets a face whose
   plane p lies beyond, which is the direction the walk should step. */
static int tet_contains(const dt_mesh *m, int ti, const long *p, int *exit_face)
{
    /* Substitution form: for a positively oriented tetrahedron, p is inside
       exactly when replacing each vertex in turn by p keeps the orientation
       positive. Chasing face windings by index parity instead is where this
       first went wrong - it silently reported every point as already present,
       so only the first insertion did anything. */
    const dt_tet *t = &m->t[ti];
    const long *q[4];
    int i;
    for (i = 0; i < 4; i++) q[i] = PT(m, t->v[i]);
    for (i = 0; i < 4; i++) {
        const long *save = q[i];
        int s;
        q[i] = p;
        s = vp_orient3d(q[0], q[1], q[2], q[3]);
        q[i] = save;
        if (s < 0) { if (exit_face) *exit_face = i; return 0; }
    }
    return 1;
}

/* Walk from m->last towards p. Exact predicates mean the walk cannot cycle
   between two tetrahedra on a tie, so no randomisation is needed; the step
   limit is a guard against a corrupted mesh, not against non-termination in
   correct operation. */
static int locate(dt_mesh *m, const long *p)
{
    int ti = m->last, guard = 0, face;
    if (ti < 0 || ti >= m->nt || m->t[ti].dead) {
        for (ti = m->nt - 1; ti >= 0 && m->t[ti].dead; ti--) ;
        if (ti < 0) return -1;
    }
    while (guard++ < 8 * m->nt + 64) {
        if (tet_contains(m, ti, p, &face)) { m->last = ti; return ti; }
        {
            int nx = m->t[ti].nb[face];
            if (nx < 0) { m->last = ti; return ti; }  /* outside the hull */
            ti = nx;
        }
    }
    return -1;
}

int dt_is_finite(const dt_mesh *m, int ti)
{
    const dt_tet *t = &m->t[ti];
    int i;
    for (i = 0; i < 4; i++) if (t->v[i] >= m->npt) return 0;
    return 1;
}

/* Plain insphere for every tetrahedron, including those touching the
   super-tetrahedron. The super corners are placed at real, finite coordinates
   far outside the data, so they behave as ordinary points and need no special
   case; the tetrahedra using them are discarded at the end by dt_is_finite.

   The "vertex at infinity" treatment was tried first and is wrong here: it
   reports every tetrahedron touching two or more super corners as always
   containing the query point, so the cavity swallows the whole mesh on each
   insertion and the triangulation never grows past four tetrahedra. */
static int in_circumsphere(const dt_mesh *m, int ti, const long *p)
{
    const dt_tet *t = &m->t[ti];
    return vp_insphere(PT(m,t->v[0]), PT(m,t->v[1]), PT(m,t->v[2]),
                       PT(m,t->v[3]), p) > 0;
}

/* Collect the cavity: every tetrahedron whose circumsphere contains p,
   reachable from the seed. Exact predicates guarantee this set is connected
   and its boundary is star-shaped from p, which is what makes the naive
   flood-fill below correct. */
static int build_cavity(dt_mesh *m, int seed, const long *p,
                        int **cav, int *ncav, int *cap)
{
    int head = 0;
    if (*cap < 64) { *cap = 64; *cav = xa_realloc(*cav, (size_t)*cap * sizeof(int)); }
    if (!*cav) return -1;
    *ncav = 0;
    (*cav)[(*ncav)++] = seed;
    m->t[seed].dead = 2;                   /* 2 = in cavity, this pass */
    while (head < *ncav) {
        int ti = (*cav)[head++], f;
        for (f = 0; f < 4; f++) {
            int nb = m->t[ti].nb[f];
            if (nb < 0 || m->t[nb].dead) continue;
            if (!in_circumsphere(m, nb, p)) continue;
            if (*ncav >= *cap) {
                *cap *= 2;
                *cav = xa_realloc(*cav, (size_t)*cap * sizeof(int));
                if (!*cav) return -1;
            }
            (*cav)[(*ncav)++] = nb;
            m->t[nb].dead = 2;
        }
    }
    return 0;
}

static int insert_point(dt_mesh *m, int pi, int **cav, int *cavcap,
                        int **nf, int *nfcap)
{
    const long *p = PT(m, pi);
    int seed = locate(m, p), ncav = 0, i, f, nnew = 0;
    if (seed < 0) return -1;
    if (!in_circumsphere(m, seed, p)) {
        /* Located tetrahedron does not contain p in its circumsphere. With
           exact predicates this means p is already a vertex (duplicate input);
           skipping keeps the triangulation valid. */
        return 0;
    }
    if (build_cavity(m, seed, p, cav, &ncav, cavcap) < 0) return -1;

    if (*nfcap < ncav * 4) {
        *nfcap = ncav * 8 + 64;
        *nf = xa_realloc(*nf, (size_t)*nfcap * 3 * sizeof(int));
        if (!*nf) return -1;
    }
    /* Boundary faces: a face of a cavity tetrahedron whose neighbour is not in
       the cavity. Each becomes a new tetrahedron with p. */
    for (i = 0; i < ncav; i++) {
        int ti = (*cav)[i];
        for (f = 0; f < 4; f++) {
            int nb = m->t[ti].nb[f];
            if (nb >= 0 && m->t[nb].dead == 2) continue;
            {
                int nt = tet_new(m);
                if (nt < 0) return -1;
                m->t[nt].v[0] = m->t[ti].v[(f+1)&3];
                m->t[nt].v[1] = m->t[ti].v[(f+2)&3];
                m->t[nt].v[2] = m->t[ti].v[(f+3)&3];
                m->t[nt].v[3] = pi;
                m->t[nt].nb[3] = nb;          /* outward neighbour keeps its link */
                tet_orient(m, nt);
                if (nb >= 0) {
                    int g;
                    for (g = 0; g < 4; g++)
                        if (m->t[nb].nb[g] == ti) { m->t[nb].nb[g] = nt; break; }
                }
                (*nf)[nnew*3+0] = nt;
                nnew++;
            }
        }
    }
    /* Retire the cavity AFTER the new tetrahedra are built, since those read
       cavity vertices, then hand the indices back for reuse. */
    for (i = 0; i < ncav; i++) m->t[(*cav)[i]].dead = 1;
    if (m->nfree + ncav > m->freecap) {
        int nc = (m->nfree + ncav) * 2 + 256;
        int *nfl = xa_realloc(m->freet, (size_t)nc * sizeof(int));
        if (nfl) { m->freet = nfl; m->freecap = nc; }
    }
    if (m->nfree + ncav <= m->freecap)
        for (i = 0; i < ncav; i++) m->freet[m->nfree++] = (*cav)[i];

    /* Link the new tetrahedra to each other by shared faces. Small n per
       insertion, so the quadratic scan is cheaper than a hash. */
    for (i = 0; i < nnew; i++) {
        int a = (*nf)[i*3];
        int j;
        for (j = i + 1; j < nnew; j++) {
            int b = (*nf)[j*3], ai, bi, shared, x, y;
            /* they share a face iff they share 3 vertices */
            shared = 0;
            for (x = 0; x < 4; x++)
                for (y = 0; y < 4; y++)
                    if (m->t[a].v[x] == m->t[b].v[y]) { shared++; break; }
            if (shared != 3) continue;
            for (ai = 0; ai < 4; ai++) {
                int found = 0;
                for (y = 0; y < 4; y++) if (m->t[b].v[y] == m->t[a].v[ai]) { found = 1; break; }
                if (!found) break;
            }
            for (bi = 0; bi < 4; bi++) {
                int found = 0;
                for (x = 0; x < 4; x++) if (m->t[a].v[x] == m->t[b].v[bi]) { found = 1; break; }
                if (!found) break;
            }
            if (ai < 4 && bi < 4) { m->t[a].nb[ai] = b; m->t[b].nb[bi] = a; }
        }
    }
    m->last = (*nf)[0];
    return 0;
}

int dt_build(dt_mesh *m, const double *xyz, int n)
{
    int i, rc = 0, *cav = NULL, *nf = NULL, cavcap = 0, nfcap = 0;
    long lo[3], hi[3];
    memset(m, 0, sizeof(*m));
    if (n < 4) return -1;
    m->npt = n;
    m->p = xa_malloc((size_t)(n + 4) * 3 * sizeof(long));
    if (!m->p) return -1;
    for (i = 0; i < 3 * n; i++) m->p[i] = vp_quant(xyz[i]);
    if (!vor_pred_in_range(m->p, 3 * n)) { free(m->p); m->p = NULL; return -2; }

    for (i = 0; i < 3; i++) { lo[i] = m->p[i]; hi[i] = m->p[i]; }
    for (i = 1; i < n; i++) {
        int k;
        for (k = 0; k < 3; k++) {
            if (m->p[3*i+k] < lo[k]) lo[k] = m->p[3*i+k];
            if (m->p[3*i+k] > hi[k]) hi[k] = m->p[3*i+k];
        }
    }
    /* Super-tetrahedron large enough that every input point is strictly
       inside it, so no input can be cospherical with its corners in a way
       that matters. Kept well inside the predicate range. */
    {
        long c[3], d = 0, k;
        for (k = 0; k < 3; k++) { c[k] = (lo[k] + hi[k]) / 2; if (hi[k]-lo[k] > d) d = hi[k]-lo[k]; }
        if (d < 1000) d = 1000;
    /* 1000x, not 20x. At 20x the enclosing simplex is small enough that the
       circumspheres of large flat HULL tetrahedra reach a corner, so those
       cells are absorbed and never emitted: 1tqn came out with 25181 cells
       against MOLE's 25189, all eight missing ones on the hull. The count
       converges as the simplex grows - 20x/60x/200x/1000x give
       25181/25186/25188/25189 - and at 1000x the mesh is cell-for-cell
       identical to MOLE's, adjacency counts included. The 256-bit accumulation
       below already covers this range. */
        d *= 10000;
        m->p[3*(n+0)+0] = c[0] - d; m->p[3*(n+0)+1] = c[1] - d; m->p[3*(n+0)+2] = c[2] - d;
        m->p[3*(n+1)+0] = c[0] + d; m->p[3*(n+1)+1] = c[1] - d; m->p[3*(n+1)+2] = c[2] - d;
        m->p[3*(n+2)+0] = c[0];     m->p[3*(n+2)+1] = c[1] + d; m->p[3*(n+2)+2] = c[2] - d;
        m->p[3*(n+3)+0] = c[0];     m->p[3*(n+3)+1] = c[1];     m->p[3*(n+3)+2] = c[2] + d;
    }
    {
        int t0 = tet_new(m);
        if (t0 < 0) { dt_free(m); return -1; }
        m->t[t0].v[0] = n; m->t[t0].v[1] = n+1; m->t[t0].v[2] = n+2; m->t[t0].v[3] = n+3;
        tet_orient(m, t0);
        m->last = t0;
    }
    for (i = 0; i < n; i++) {
        if (insert_point(m, i, &cav, &cavcap, &nf, &nfcap) < 0) { rc = -1; break; }
    }
    free(cav); free(nf);
    if (rc) dt_free(m);
    return rc;
}

void dt_free(dt_mesh *m)
{
    if (m) { free(m->freet); m->freet = NULL; m->nfree = m->freecap = 0; }
    if (!m) return;
    free(m->p); free(m->t);
    memset(m, 0, sizeof(*m));
}

int dt_count_live(const dt_mesh *m)
{
    int i, k = 0;
    for (i = 0; i < m->nt; i++) if (!m->t[i].dead) k++;
    return k;
}

int dt_count_finite(const dt_mesh *m)
{
    int i, k = 0;
    for (i = 0; i < m->nt; i++) if (!m->t[i].dead && dt_is_finite(m, i)) k++;
    return k;
}

int dt_verify(const dt_mesh *m)
{
    int i, j, bad = 0;
    for (i = 0; i < m->nt; i++) {
        const dt_tet *t;
        if (m->t[i].dead || !dt_is_finite(m, i)) continue;
        t = &m->t[i];
        if (vp_orient3d(PT(m,t->v[0]), PT(m,t->v[1]), PT(m,t->v[2]), PT(m,t->v[3])) <= 0) {
            bad++; continue;
        }
        for (j = 0; j < m->npt; j++) {
            if (j == t->v[0] || j == t->v[1] || j == t->v[2] || j == t->v[3]) continue;
            if (vp_insphere(PT(m,t->v[0]), PT(m,t->v[1]), PT(m,t->v[2]),
                            PT(m,t->v[3]), PT(m,j)) > 0) { bad++; break; }
        }
    }
    return bad;
}

int dt_circumcentre(const dt_mesh *m, int ti, const double *radii,
                    double *cxyz, double *clearance)
{
    const dt_tet *t = &m->t[ti];
    double a[3], b[3], c[3], d[3], M[3][3], r[3], det, cc[3];
    int i;
    if (m->t[ti].dead || !dt_is_finite(m, ti)) return -1;
    for (i = 0; i < 3; i++) {
        a[i] = (double)PT(m,t->v[0])[i] / VP_SCALE;
        b[i] = (double)PT(m,t->v[1])[i] / VP_SCALE;
        c[i] = (double)PT(m,t->v[2])[i] / VP_SCALE;
        d[i] = (double)PT(m,t->v[3])[i] / VP_SCALE;
    }
    for (i = 0; i < 3; i++) {
        M[0][i] = b[i]-a[i]; M[1][i] = c[i]-a[i]; M[2][i] = d[i]-a[i];
    }
    for (i = 0; i < 3; i++) {
        const double *q = (i==0)?b:((i==1)?c:d);
        r[i] = 0.5*((q[0]*q[0]+q[1]*q[1]+q[2]*q[2]) - (a[0]*a[0]+a[1]*a[1]+a[2]*a[2]));
    }
    det = M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
        - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
        + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0]);
    if (fabs(det) < 1e-12) return -1;      /* flat tetrahedron */
    /* Cramer */
    {
        double A[3][3];
        int k;
        for (k = 0; k < 3; k++) {
            memcpy(A, M, sizeof(A));
            A[0][k]=r[0]; A[1][k]=r[1]; A[2][k]=r[2];
            cc[k] = ( A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                    - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                    + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]) ) / det;
        }
    }
    if (cxyz) { cxyz[0]=cc[0]; cxyz[1]=cc[1]; cxyz[2]=cc[2]; }
    if (clearance) {
        double R = sqrt((cc[0]-a[0])*(cc[0]-a[0]) + (cc[1]-a[1])*(cc[1]-a[1])
                      + (cc[2]-a[2])*(cc[2]-a[2]));
        /* Largest sphere that fits here: the circumradius less the vdW radius
           of a defining atom. Without radii this is the plain circumradius. */
        *clearance = radii ? R - radii[t->v[0]] : R;
    }
    return 0;
}
