/* MOLE 2's own Delaunay triangulation (DHTriangulation). See mole_dh.c. */
#ifndef MOLE_DH_H
#define MOLE_DH_H

/* Finite cells only, in MOLE's own cell order, with MOLE's own per-cell vertex
   order. tv[4t+k] indexes the input points; tn[4t+k] is the neighbouring cell
   across the face OPPOSITE tv[4t+k], or -1 when that neighbour is infinite. */
typedef struct {
    int  nt;
    int *tv, *tn;
} dh_mesh;

/* Returns 0 on success. n must be >= 4. */
int  dh_build(dh_mesh *out, const double *xyz, int n);
void dh_mesh_free(dh_mesh *m);

/* Hilbert index of an integer cell, exposed for the oracle test. */
int  dh_hilbert_encode(int order, int x, int y, int z);

#endif
