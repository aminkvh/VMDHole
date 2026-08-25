/* MOLE 2's cavity pipeline. See mole_complex.c. */
#ifndef MOLE_COMPLEX_H
#define MOLE_COMPLEX_H

#include <stdio.h>
#include "../voronoi/vor_delaunay.h"

/* TunnelWeightFunction. VoronoiScale is MOLE's default (enum value 0, and set
   explicitly by ComplexParameters); LengthAndRadius is the classic
   int ds/rho^2 that CAVER and MOLE 1.x use. */
typedef enum {
    MOLE_W_VORONOI_SCALE = 0,
    MOLE_W_LENGTH_AND_RADIUS,
    MOLE_W_LENGTH,
    MOLE_W_CONSTANT
} mole_weight_fn;

typedef struct {
    double probe_radius;        /* MOLE ProbeRadius,       default 3    */
    double interior_threshold;  /* MOLE InteriorThreshold, default 1.25 */
    int    min_depth;           /* MOLE MinDepth,          default 8    */
    double min_depth_length;    /* MOLE MinDepthLength,    default 5    */
    double min_tunnel_length;   /* MOLE MinTunnelLength,   default 0    */
    mole_weight_fn weight;      /* MOLE WeightFunction,    default VoronoiScale */
    /* MOLE StrictInterior, default false. Degenerate upstream and therefore
       here: it makes MaxClearance zero everywhere, so the interior removal
       takes every tetrahedron and no cavity, and hence no tunnel, survives.
       See mole_complex.c. */
    int    strict_interior;
} mole_params;

#define MOLE_MAXATOM 200000

/* One structure's atoms, in file order. */
typedef struct {
    int    n;
    double xyz[3*MOLE_MAXATOM];
    char   elem[MOLE_MAXATOM][8], chain[MOLE_MAXATOM][8], resn[MOLE_MAXATOM][8];
    int    water[MOLE_MAXATOM], seq[MOLE_MAXATOM];
    /* Optional trailing columns: B-factor (BRadius), then the PDB atom name.
       backbone and freeatom are DERIVED from the name by MOLE's own rules
       rather than supplied, so the caller cannot disagree with MOLE about what
       a backbone atom is - two callers already had. has_names records whether
       the name column was there; without it FreeRadius has no atoms to measure
       and is reported as Radius, and every layer counts as non-backbone. */
    double bfac[MOLE_MAXATOM];
    char   name[MOLE_MAXATOM][8];
    int    backbone[MOLE_MAXATOM], freeatom[MOLE_MAXATOM], has_names;
    /* Columns 11 and 12: alternate-location indicator and residue insertion
       code, "-" when blank. They do not change the file order - which is what
       the jitter is drawn in - but they DO change the PIVOT order, and the
       pivot order is the triangulation's insertion history. See mole_pivots. */
    char   alt[MOLE_MAXATOM][4], icode[MOLE_MAXATOM][4];
    /* Pivot order, filled by mole_pivots: order[p] is the atom index of pivot p.
       mole_pivot_extras and mole_pivot_residues MUST walk this, not file order,
       or their per-pivot arrays index a different atom than piv[] does. */
    int    order[MOLE_MAXATOM], npiv;
    /* Pivot coordinates BEFORE the general-position jitter, in pivot order.
       FindHetResidues searches Structure.InvariantKdAtomTree(), and
       InvariantPosition is set at construction and never touched by the jitter
       - so that one search runs on different numbers from everything else. */
    double invpiv[3*MOLE_MAXATOM];
} mole_atoms;

typedef struct {
    int    nt;                  /* finite tetrahedra */
    const double *axyz;         /* atom positions, jittered */
    const double *arad;         /* atom vdW radii */
    int    *tv, *tn;            /* 4 vertices / 4 neighbours per tetrahedron */
    double *center, *vcenter;   /* centroid and circumcentre */
    double *volume, *maxclear;
    int    *depth;
    double *depthlen;
    double *eclear, *elen, *eweight;   /* per (tetrahedron, face) */
    double *evweight;                  /* VoronoiScale weight - MOLE's DEFAULT */
    char   *boundary, *alive;
    /* SurfaceCavity membership: `alive` as it stands at the snapshot between the
       probe peel and the interior removal. Retained because custom exits and
       paths search that graph, and it is not any cavity's - MOLE keeps it as a
       cavity in its own right with Id 0. */
    char   *surface;
    int    *comp;
    /* cavityGraph.Vertices order: the position at which each tetrahedron first
       appears while AddVerticesAndEdgeRange walks Triangulation.Edges. MOLE
       filters THIS order to get a cavity's boundary tetrahedra, and that order
       decides which opening pivots come out. See mole_openings. */
    int    *vorder;
    int     n_surface;          /* size of the SurfaceCavity snapshot */
} mole_complex;

typedef struct {
    int    count, depth, has_boundary;
    double volume, depth_length;
} mole_cavity;

double mole_vdw_radius(const char *elem);
/* PdbEx.IsBackboneAtom: an amino OR nucleic backbone atom NAME. */
int    mole_is_backbone_name(const char *name);
/* PdbResidue.IsAminoName: one of the 20 standard residues. */
int    mole_is_amino_name(const char *resn);
/* Chain-identifier ORDER, for every OrderBy(ChainIdentifier) MOLE performs.
   InvariantCulture collation, not strcmp - see mole_complex.c. Equality tests
   stay on strcmp; only ordering goes through here. */
int    mole_chain_cmp(const char *a, const char *b);
/* MOLE's decimal parser, reproduced ULP-for-ULP. See mole_complex.c. */
double mole_parse_double(const char *s);
double mole_det4(const double m[4][4]);
/* MathHelper.SphereFromPoints' centre. Exposed so the oracle test can pin it to
   MOLE's own measured circumcentre rather than to our compiled result. */
void   mole_circumsphere(const double p[4][3], double *cx, double *cy, double *cz);
int    mole_read_atoms(const char *path, mole_atoms *A, int cap);
int    mole_pivots(mole_atoms *A, double *piv, double *rad, int *pivsrc);
/* Per-pivot B-factor, free flag and backbone flag, in pivot order. Any output
   may be NULL. Call after mole_pivots. */
int    mole_pivot_extras(const mole_atoms *A, double *bfac, int *freeatom,
                         int *backbone);
int    mole_build_from_tetra(mole_complex *M, const int *tv, const int *tn, int nt,
                             const double *axyz, const double *arad,
                             const mole_params *P);
int    mole_build(mole_complex *M, const dt_mesh *m, const double *axyz,
                  const double *arad, const mole_params *P);
int    mole_cavities(mole_complex *M, const mole_params *P,
                     mole_cavity **out, int *n_channel, int *n_void);
/* PdbEx.IsHetAtom, as far as an atom table can express it: MOLE's rule is
   "HETATM record OR not a standard residue name", and a table has no record
   type, so only the second half is testable. */
int    mole_is_het_resname(const char *resn);
void   mole_edge_update(mole_complex *M, int t, int k);
void   mole_free(mole_complex *M);

#endif
