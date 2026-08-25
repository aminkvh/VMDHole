/* Exact geometric predicates for the Voronoi/Delaunay tunnel search.
 *
 * WHY EXACT, AND WHY NOT FLOATING POINT
 * -------------------------------------
 * Incremental Delaunay does not use the VALUE of orient3d/insphere, only its
 * SIGN. A sign that is wrong near a degeneracy does not perturb the result
 * slightly - it corrupts the topology: point-location walks that never
 * terminate, tetrahedra with inverted orientation, a mesh that is not a
 * triangulation of anything. Protein coordinates produce near-degenerate
 * configurations readily (crystallographic symmetry, repeated geometry,
 * coordinates written to a fixed grid), so this is not a corner case here.
 *
 * WHY INTEGERS RATHER THAN SHEWCHUK EXPANSIONS
 * --------------------------------------------
 * The usual remedy is adaptive floating-point expansion arithmetic. It is a
 * lot of delicate code whose failure mode is a silently wrong sign. PDB
 * coordinates carry exactly three decimals, so scaling by 1000 is exact and
 * turns both predicates into integer determinants, which __int128 evaluates
 * with no rounding at all. Degenerate inputs then return exactly 0 by
 * construction rather than "small".
 *
 * Quantisation to 0.001 A is a real but harmless approximation: it perturbs a
 * coordinate by at most 0.0005 A, far below any radius this program reports,
 * and it buys exactness and reproducibility. Trajectory coordinates that carry
 * more precision than a PDB are quantised the same way, which also makes the
 * triangulation reproducible across runs.
 *
 * RANGE
 * -----
 * With SCALE = 1000 and |coord| <= 1e6 scaled (1000 A), the largest quantity
 * is insphere's degree-5 term: (2e6)^2 * (2e6)^3 = 3.2e31, against __int128's
 * ~1.7e38. orient3d peaks at (2e6)^3 = 8e18. Both comfortable.
 * vor_pred_in_range() rejects anything outside that envelope rather than
 * overflowing silently.
 *
 * WHAT THE TESTS DO AND DO NOT SHOW
 * ---------------------------------
 * vor_predicates_test.c confirms exactness directly: orient3d reproduces a
 * closed-form determinant (-3t) on 20000 near-degenerate inputs, and both
 * predicates return exactly 0 on constructed coplanar/cospherical sets.
 *
 * It does NOT show naive double failing at this coordinate scale, and the
 * attempt is recorded rather than dropped. orient3d's 2x2 minors stay under
 * 2^53 for |coord| <= 1e6, so double evaluates them exactly; the axis-aligned
 * cospherical cases likewise contain enough zeros to stay exact. So the case
 * for exactness here is not "double is observably wrong on protein data" - it
 * is that the failure mode is a wrong SIGN producing corrupt topology rather
 * than a slightly wrong number, and that it is cheap to remove entirely.
 */
#ifndef VOR_PREDICATES_H
#define VOR_PREDICATES_H

typedef __int128 vp_int;

/* Overridable at build time so a second triangulator can be compiled from these
   same sources at a finer grid. The MOLE port needs one: MOLE jitters atoms by
   +-0.00005 A to force general position, which at SCALE 1000 is +-0.1 of a
   quantisation step and is rounded away entirely, putting the atoms back on the
   degenerate grid the jitter exists to escape. SCALE 1e5 keeps the jitter at
   +-5 steps; insphere's degree-5 term then peaks at 2.7e35 for a 61 A protein
   against __int128's 1.7e38. SCALE 1e6 overflows, so 1e5 is the only workable
   setting and it is not a free choice. Default is unchanged, so every existing
   caller keeps its current output. */
#ifndef VP_SCALE
#define VP_SCALE      1000.0
#endif
#ifndef VP_MAX_COORD
#define VP_MAX_COORD  1000000L   /* 1000 A, scaled */
#endif

/* Quantise one coordinate to the exact integer grid. Round-half-away-from-zero
   so the mapping is symmetric about 0 and independent of FP rounding mode. */
long vp_quant(double x);

/* 1 if every coordinate is representable without overflow, 0 otherwise. */
int vor_pred_in_range(const long *p, int n3);

/* Sign of the orientation determinant of a,b,c,d (each 3 scaled longs):
 *   > 0  d lies BELOW the plane abc  (abc counterclockwise seen from d)
 *   = 0  the four points are exactly coplanar
 *   < 0  d lies ABOVE
 * Exact: the returned sign is the true sign of the determinant. */
int vp_orient3d(const long *a, const long *b, const long *c, const long *d);

/* Sign of the insphere determinant of a,b,c,d,e (each 3 scaled longs):
 *   > 0  e lies INSIDE the sphere through a,b,c,d
 *   = 0  the five points are exactly cospherical
 *   < 0  e lies outside
 * Assumes orient3d(a,b,c,d) > 0; the caller must orient the tetrahedron first.
 * Exact. */
int vp_insphere(const long *a, const long *b, const long *c,
                const long *d, const long *e);

#endif
