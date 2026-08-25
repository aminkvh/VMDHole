/* MOLE 2's lining layers and physico-chemical properties. See mole_lining.c. */
#ifndef MOLE_LINING_H
#define MOLE_LINING_H

#include "mole_tunnel.h"

#define MOLE_MAXRES 60000
#define MOLE_LAYER_MAXRES 5      /* at most five atoms, so at most five residues */

/* Unique residues among the pivots, in first-seen (file) order. */
typedef struct {
    int  n;
    char name[MOLE_MAXRES][8], chain[MOLE_MAXRES][8];
    int  seq[MOLE_MAXRES];
} mole_residues;

/* TunnelPhysicoChemicalProperties. Mutability is an int because MOLE casts it. */
typedef struct {
    int    charge, ionizable, npos, nneg, mutability;
    double hydropathy, hydrophobicity, polarity, logp, logd, logs;
} mole_props;

/* TunnelLayer, after the equal-lining runs have been merged. */
/* TunnelLining.ResidueFlow: every (residue, backbone-touched) pair the tunnel
   meets, in the order it first meets them. The SAME residue appears twice when
   it lines one layer at its backbone and another at its side chain. */
typedef struct { int res; char isbb; } mole_flow_entry;

typedef struct {
    int    nres;
    int    res[MOLE_LAYER_MAXRES];    /* index into mole_residues */
    char   isbb[MOLE_LAYER_MAXRES];   /* residue touched ONLY at backbone */
    int    flow[MOLE_LAYER_MAXRES];   /* index into mole_lining.flow */
    double start, end;                /* Start/EndDistance */
    double radius, freeradius, bradius;   /* MIN over the merged run */
    int    localmin;
    mole_props props;
} mole_layer;

typedef struct {
    int              nl;
    mole_layer      *layer;
    int              nflow;
    mole_flow_entry *flow;
    mole_props       props;   /* whole tunnel, unweighted */
    mole_props       wprops;  /* weighted by layer length - what tunnels.csv prints */
} mole_lining;

/* Per-pivot residue index plus the residue table. Call after mole_pivots. */
int  mole_pivot_residues(const mole_atoms *A, int *pres, mole_residues *R);
/* Build the layers of one tunnel. pbb/pres are per pivot, R the residue table. */
int  mole_lining_build(const mole_tunnel_profile *p,
                       const double *axyz, const double *arad, int na,
                       const int *pres, const mole_residues *R, const int *pbb,
                       mole_lining *out);
void mole_lining_free(mole_lining *L);
/* TunnelLining.ComputeFlow, exposed for direct testing of the chain-ordering
   regression (mole_lining_regression_test.c). See mole_lining.c. */
int  mole_build_flow(mole_layer *L, int nl, const mole_residues *R,
                    mole_flow_entry *flow);
/* Cavity.Create's boundary/inner residue split, both sorted by chain then
   number (the residue table is already in file order, which is that order for a
   single-chain structure; the caller sorts). Returns (nboundary<<16)|ninner. */
int  mole_cavity_residues(const mole_complex *M, int comp, const int *pres,
                          const mole_residues *R, int *bnd, int *inner);
/* CalculateResidueProperties: the RESIDUES overload, where every residue also
   contributes the backbone constants. NOT the per-layer one. */
void mole_cavity_properties(const int *res, int nres, const mole_residues *R,
                            mole_props *P);

/* Tunnel.FindHetResidues - the HET residues the tunnel passes through, as
   <HetResidues> reports them. Returns the count; `out` holds residue indices
   into R, ordered by chain then number. See mole_lining.c. */
int mole_het_residues(const mole_tunnel_profile *pr, const mole_atoms *A,
                      const double *piv, int np,
                      const mole_complex *M, const int *path, int npath,
                      const int *pres, const mole_residues *R, int *out, int cap);
/* The OrderBy(ChainIdentifier).ThenBy(Number) tail of mole_het_residues,
   exposed for direct testing with a synthetic residue table. */
void mole_sort_residues_by_chain(int *out, int nout, const mole_residues *R);

#endif
