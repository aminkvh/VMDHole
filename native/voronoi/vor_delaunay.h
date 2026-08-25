/* Incremental 3-D Delaunay triangulation, built on the exact predicates in
 * vor_predicates.h.
 *
 * WHAT THIS IS FOR
 * ----------------
 * The Voronoi vertices of a set of atom centres are the circumcentres of the
 * Delaunay tetrahedra, and the Voronoi edges connect circumcentres of
 * face-adjacent tetrahedra. That edge network is the search graph CAVER uses,
 * and it is what the lattice search approximates. Building the Delaunay
 * triangulation is therefore the whole of the geometric work; the tunnel search
 * on top of it is a shortest-path problem on the resulting graph.
 *
 * DESIGN NOTES
 * ------------
 * Bowyer-Watson insertion with a walking point location. Every geometric
 * decision goes through vp_orient3d / vp_insphere, which are exact, so the
 * usual failure modes of a floating-point implementation - a cavity that is not
 * star-shaped, a walk that cycles, tetrahedra with inverted orientation - are
 * excluded by construction rather than guarded against with tolerances.
 *
 * Tetrahedra are stored with their four vertices and their four neighbours,
 * where neighbour[i] is the tetrahedron opposite vertex[i]. Deleted tetrahedra
 * are marked rather than compacted, so indices stay stable while a cavity is
 * being retriangulated.
 *
 * The four vertices of the enclosing super-tetrahedron are indices
 * -1..-4 (stored as n..n+3 internally). A tetrahedron touching one of them is
 * outside the convex hull; callers that want only interior Voronoi vertices
 * must skip those, and dt_is_finite() answers that.
 */
#ifndef VOR_DELAUNAY_H
#define VOR_DELAUNAY_H

#include "vor_predicates.h"

typedef struct {
    int v[4];        /* vertex indices; >= npt means a super-tetrahedron corner */
    int nb[4];       /* nb[i] = tetra opposite v[i], or -1 */
    int dead;        /* 1 once removed by an insertion */
} dt_tet;

typedef struct {
    long  *p;        /* 3 * (npt + 4) quantised coordinates */
    int    npt;      /* real points; the 4 super corners follow */
    dt_tet *t;
    int    nt, cap;
    int    last;     /* where the previous walk ended, as the next walk's start */
    /* Tetrahedra retired by an insertion, for reuse. Without this the array
       only grows: a 12-ball build reached 1.88M entries of which 1.6M were
       dead. Reclaiming them holds it at 274k. NOTE this is a memory fix, not a
       speed one - it was measured to change build time not at all, because the
       bottleneck is seg_clearance over edges, not insertion. */
    int   *freet;
    int    nfree, freecap;
} dt_mesh;

/* Build from n points given as 3n doubles (Angstroms). Returns 0 on success.
   Coordinates are quantised by vp_quant; out-of-range input is rejected. */
int  dt_build(dt_mesh *m, const double *xyz, int n);
void dt_free(dt_mesh *m);

/* 1 if the tetrahedron uses only real points (no super-tetrahedron corner). */
int  dt_is_finite(const dt_mesh *m, int ti);

/* Number of live tetrahedra, and of live finite ones. */
int  dt_count_live(const dt_mesh *m);
int  dt_count_finite(const dt_mesh *m);

/* Verification, for tests and for a paranoid caller:
   returns 0 if every live finite tetrahedron is positively oriented AND no
   input point lies strictly inside any live finite tetrahedron's circumsphere.
   Non-zero is the number of violations found. O(nt * npt), so it is a test
   tool, not something to run in production. */
int  dt_verify(const dt_mesh *m);

/* Circumcentre of tetrahedron ti in Angstroms, and the radius of the largest
   empty sphere there (circumradius minus the vdW radius of a defining atom,
   if radii are supplied; pass NULL for plain circumradius). Returns 0 on
   success, non-zero if the tetrahedron is degenerate. */
int  dt_circumcentre(const dt_mesh *m, int ti, const double *radii,
                     double *cxyz, double *clearance);

#endif
