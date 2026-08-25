/* MOLE 2's tunnel stage. See mole_tunnel.c. */
#ifndef MOLE_TUNNEL_H
#define MOLE_TUNNEL_H

#include "mole_complex.h"

typedef struct { int n; double *t, *y, *m; } mole_spline;

typedef struct {
    mole_spline sx, sy, sz, sr;
    /* MOLE's two extra widths, on the same 100 samples as sr:
       FreeRadius ignores side chains (measured against backbone + het only),
       BRadius adds the mean B-factor RMSF of the same five nearest atoms. */
    mole_spline sfr, sbr;
    double length;
} mole_tunnel_profile;

int  mole_openings(const mole_complex *M, int cavity_comp, double cover_radius,
                   int **out);
/* Cavity.GetOpening - snap a user point to a cavity's nearest boundary-facet
   tetrahedron by CENTROID, rejected past radius. mask selects the cavity
   (surface membership, or NULL for the ordinary alive/comp test). */
int  mole_cavity_opening(const mole_complex *M, int cavity_comp, const char *mask,
                         const double *point, double radius);
/* Capacity of the caller's `out` array, and therefore the hard upper bound on
   max_origins. The caller's buffer is a fixed stack array, so an unclamped
   user value (the GUI accepts any nonnegative integer) plus a cavity with
   enough qualifying maxima would write past it. Both the caller and
   mole_auto_origins itself clamp to this - belt and braces, because the
   overflow is silent and corrupts the stack. */
#define MOLE_MAX_ORIGINS 16
int  mole_auto_origins(const mole_complex *M, int cavity_comp,
                       double cover_radius, int max_origins, int *out);
double mole_edge_cost(const mole_complex *M, int t, int k, mole_weight_fn w);
int  mole_dijkstra(const mole_complex *M, int cavity_comp, int src,
                   double *dist, int *prev, mole_weight_fn w);
/* The same, but a non-NULL mask REPLACES the alive/comp membership test - used
   to search the SurfaceCavity, whose tetrahedra are neither alive nor in any
   cavity component by then. */
int  mole_dijkstra_mask(const mole_complex *M, int cavity_comp, int src,
                        double *dist, int *prev, mole_weight_fn w,
                        const char *mask);
void mole_spline_init(mole_spline *s, const double *t, const double *y, int n);
double mole_spline_eval(const mole_spline *s, double x);
void mole_spline_free(mole_spline *s);
double mole_radius_at_raw(const double *p, const double *axyz, const double *arad, int na);
/* The five nearest atoms to p, nearest first; returns how many (<5 only when
   the structure has fewer atoms). sel must hold 5. */
int  mole_nearest5(const double *p, const double *axyz, const double *arad,
                   int na, int *sel);
int  mole_control_path(const mole_complex *M, const int *path, int np,
                       const double *axyz, const double *arad, int na,
                       double interior_threshold, int *out);
/* is_path = TunnelType.Path: keeps the first tetrahedron regardless of radius,
   which is what CalculateProfile's `Type == Tunnel` guard means. */
int  mole_control_path_ex(const mole_complex *M, const int *path, int np,
                          const double *axyz, const double *arad, int na,
                          double interior_threshold, int *out, int is_path);
/* Cavity.GetTetrahedron - nearest member by CENTROID over the WHOLE cavity,
   not just its boundary facets. Used to snap a path's two endpoints. */
int  mole_cavity_tetrahedron(const mole_complex *M, int cavity_comp,
                             const char *mask, const double *point, double radius);
int  mole_profile(const mole_complex *M, const int *path, int np,
                  const double *axyz, const double *arad, int na,
                  mole_tunnel_profile *out);
/* FreeRadius/BRadius need per-atom B-factors and the backbone-or-het mask.
   Pass NULL for either to leave that width equal to the plain radius. */
void mole_profile_extras(const double *bfac, const int *freeatom);
void mole_profile_free(mole_tunnel_profile *p);
void mole_filter_similar(const mole_tunnel_profile *prof, const int *path_len,
                         char *dead, int n, double max_similarity);
int  mole_filter_bottleneck(const mole_tunnel_profile *p, double bottleneck,
                            double tolerance, double density);

#endif
