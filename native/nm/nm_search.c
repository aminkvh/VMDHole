/* nm_search.c - prototype: HOLE-style pore profile via Nelder-Mead.
 *
 * Goal-1 feasibility probe for the <30ms/frame target: replace HOLE's
 * per-slice simulated annealing (1000 MC steps) with a Nelder-Mead simplex
 * seeded from the previous slice (CHAP's published approach), on top of a
 * topology-loaded-once architecture (radii assigned one time, coordinates
 * streamable).
 *
 *   cc -O2 -o nm_search nm_search.c -lm
 *   ./nm_search frame.pdb hole.rad  cx cy cz  vx vy vz  sample endrad
 *
 * Prints one line per slice:  t x y z r   (t = coordinate along cvect)
 * and timing to stderr.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nm_threads.h"
#include "../hole_io.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXATOM 200000

static double ax[MAXATOM], ay[MAXATOM], az[MAXATOM], ar[MAXATOM];
static char aresn[MAXATOM][4];   /* residue name, for --ignore */
static int natom = 0;

/* ---- radius table ---------------------------------------------------- */
#define MAXRAD 64
static char rn[MAXRAD][5], rres[MAXRAD][4];
static double rval[MAXRAD];
static int nrad = 0;

static int wmatch(const char *pat, const char *s, int n) {
    for (int i = 0; i < n; i++)
        if (pat[i] != '?' && pat[i] != s[i]) return 0;
    return 1;
}

static void load_rad(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    char line[256];
    while (fgets(line, sizeof line, f)) {
        char a[16], b[16]; double v;
        if (sscanf(line, "VDWR %15s %15s %lf", a, b, &v) == 3 && nrad < MAXRAD) {
            snprintf(rn[nrad], 5, "%-4.4s", a);
            snprintf(rres[nrad], 4, "%-3.3s", b);
            rval[nrad++] = v;
        }
    }
    fclose(f);
}

static double radius_for(const char *name4, const char *res3) {
    for (int i = 0; i < nrad; i++)
        if (wmatch(rn[i], name4, 4) && wmatch(rres[i], res3, 3))
            return rval[i];
    return -1.0;
}

/* ---- pdb ------------------------------------------------------------- */
/* The plugin's packed coordinate record (input_frame.vhb, written by
   _write_hole_coord_bin): "VMDHOLEC", int32 version, int32 natoms, then
   5*n name field (PDB cols 12-16), 3*n residue names, n chain ids, int32*n
   resids, and x[n], y[n], z[n] doubles, all native-endian. Coordinates are
   rounded through "%.3f" so this path gives exactly the result the PDB path
   gives for the same frame. */
static int load_vhb(FILE *f) {
    int hdr[2];
    if (fread(hdr, sizeof(int), 2, f) != 2 || hdr[0] != 1 || hdr[1] < 1 || hdr[1] >= MAXATOM) return 0;
    int n = hdr[1];
    char *n5 = malloc(5 * (size_t)n), *r3 = malloc(3 * (size_t)n), *c1 = malloc((size_t)n);
    int *rid = malloc(sizeof(int) * (size_t)n);
    double *xyz = malloc(3 * sizeof(double) * (size_t)n);
    int ok = fread(n5, 5, n, f) == (size_t)n && fread(r3, 3, n, f) == (size_t)n &&
             fread(c1, 1, n, f) == (size_t)n && fread(rid, sizeof(int), n, f) == (size_t)n &&
             fread(xyz, sizeof(double), 3 * (size_t)n, f) == 3 * (size_t)n;
    if (ok) {
        for (int i = 0; i < n; i++) {
            char name[5] = "    ", res[4] = "   ";
            int k = 0;
            for (int j = 0; j < 4; j++) if (n5[5*i+j] != ' ') name[k++] = n5[5*i+j];
            memcpy(res, r3 + 3*i, 3);
            double r = radius_for(name, res);
            if (r < 0) { fprintf(stderr, "no radius: '%s' '%s'\n", name, res); exit(1); }
            char buf[32];
            snprintf(buf, sizeof buf, "%.3f", xyz[i]);       ax[i] = atof(buf);
            snprintf(buf, sizeof buf, "%.3f", xyz[n + i]);   ay[i] = atof(buf);
            snprintf(buf, sizeof buf, "%.3f", xyz[2*n + i]); az[i] = atof(buf);
            ar[i] = r;
            memcpy(aresn[i], res, 3); aresn[i][3] = 0;
        }
        natom = n;
    }
    free(n5); free(r3); free(c1); free(rid); free(xyz);
    return ok;
}

static void load_pdb(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(1); }
    char line[256];
    if (fread(line, 1, 8, f) == 8 && !memcmp(line, "VMDHOLEC", 8)) {
        if (!load_vhb(f)) { fprintf(stderr, "%s: unreadable packed coordinate record\n", path); exit(1); }
        fclose(f);
        return;
    }
    fclose(f);
    hio_reader rd; hio_rec a;
    if (!hio_open(&rd, path)) { perror(path); exit(1); }
    while (hio_next(&rd, &a)) {
        char name[5] = "    ";
        /* HOLE left-justifies the stripped atom name before matching */
        int k = 0;
        for (int i = 0; i < 4; i++)
            if (a.name[i] != ' ') name[k++] = a.name[i];
        double r = radius_for(name, a.resn);
        if (r < 0) { fprintf(stderr, "no radius: '%s' '%s'\n", name, a.resn); exit(1); }
        ax[natom] = a.x; ay[natom] = a.y; az[natom] = a.z;
        ar[natom] = r;
        memcpy(aresn[natom], a.resn, 4);
        if (++natom >= MAXATOM) { fprintf(stderr, "too many atoms\n"); exit(1); }
    }
    hio_close(&rd);
}

/* streaming: re-read COORDINATES only (topology/radii keep the first
 * frame's assignment; atom order must match, which VMD's writepdb of a
 * fixed selection guarantees). Returns 0 on atom-count mismatch. */
/* ---- cell grid over atoms -------------------------------------------- */
static double cell = 4.0, gx0, gy0, gz0, maxr = 0;
static int nx, ny, nz;
static int *cellstart, *cellatom;
/* Atom x,y,z,r laid out in CELL ORDER, four doubles together. The clearance
   loop reads all four of an atom and then moves to the next atom in the same
   cell, so cell-ordered and interleaved is one sequential stream instead of
   four scattered lookups through cellatom. */
static double *cellxyzr;

static void build_grid(void) {
    /* Cell size is the one acceleration parameter here: too small and the
       expanding-ring search walks many empty cells to reach a wide pore's
       nearest atom, too large and each cell holds atoms the ring bound can no
       longer exclude. Overridable so it can be swept per system rather than
       assumed; the result is a pure speed choice, never a different answer. */
    { const char *e = getenv("NM_CELL");
      if (e && atof(e) > 0.5) cell = atof(e); }
    free(cellstart); free(cellatom); free(cellxyzr);   /* streaming rebuilds per frame */
    double xlo = 1e30, ylo = 1e30, zlo = 1e30, xhi = -1e30, yhi = -1e30, zhi = -1e30;
    maxr = 0;
    for (int i = 0; i < natom; i++) {
        if (ax[i] < xlo) xlo = ax[i]; if (ax[i] > xhi) xhi = ax[i];
        if (ay[i] < ylo) ylo = ay[i]; if (ay[i] > yhi) yhi = ay[i];
        if (az[i] < zlo) zlo = az[i]; if (az[i] > zhi) zhi = az[i];
        if (ar[i] > maxr) maxr = ar[i];
    }
    gx0 = xlo - 1e-6; gy0 = ylo - 1e-6; gz0 = zlo - 1e-6;
    nx = (int)((xhi - gx0) / cell) + 1;
    ny = (int)((yhi - gy0) / cell) + 1;
    nz = (int)((zhi - gz0) / cell) + 1;
    int ncell = nx * ny * nz;
    cellstart = calloc(ncell + 1, sizeof(int));
    cellatom = malloc(natom * sizeof(int));
    for (int i = 0; i < natom; i++) {
        int c = ((int)((ax[i] - gx0) / cell)) * ny * nz
              + ((int)((ay[i] - gy0) / cell)) * nz
              + (int)((az[i] - gz0) / cell);
        cellstart[c + 1]++;
    }
    for (int c = 0; c < ncell; c++) cellstart[c + 1] += cellstart[c];
    int *fill = calloc(ncell, sizeof(int));
    for (int i = 0; i < natom; i++) {
        int c = ((int)((ax[i] - gx0) / cell)) * ny * nz
              + ((int)((ay[i] - gy0) / cell)) * nz
              + (int)((az[i] - gz0) / cell);
        cellatom[cellstart[c] + fill[c]++] = i;
    }
    free(fill);
    cellxyzr = malloc((size_t)natom * 4 * sizeof(double));
    for (int k = 0; k < natom; k++) {
        int i = cellatom[k];
        cellxyzr[4*k+0] = ax[i]; cellxyzr[4*k+1] = ay[i];
        cellxyzr[4*k+2] = az[i]; cellxyzr[4*k+3] = ar[i];
    }
}

/* clearance = min over atoms of (dist - vdw radius); expanding-ring search
 * with early exit once the ring's lower bound exceeds the current best */
static double clearance(double px, double py, double pz) {
    int cx = (int)((px - gx0) / cell), cy = (int)((py - gy0) / cell),
        cz = (int)((pz - gz0) / cell);
    double best = 1e30;
    /* enough rings to reach every cell even when p is outside the grid */
    int mx = (cx < 0 ? nx - cx : (cx >= nx ? cx : (cx > nx - 1 - cx ? cx : nx - 1 - cx)));
    int my = (cy < 0 ? ny - cy : (cy >= ny ? cy : (cy > ny - 1 - cy ? cy : ny - 1 - cy)));
    int mz = (cz < 0 ? nz - cz : (cz >= nz ? cz : (cz > nz - 1 - cz ? cz : nz - 1 - cz)));
    int maxring = (mx > my ? (mx > mz ? mx : mz) : (my > mz ? my : mz));
    for (int ring = 0; ring <= maxring; ring++) {
        double lb = (ring - 1) * cell - maxr;   /* closest any atom in this ring can be */
        if (ring > 0 && lb > best) break;
        int xlo = cx - ring, xhi = cx + ring;
        for (int ix = xlo; ix <= xhi; ix++) {
            if (ix < 0 || ix >= nx) continue;
            for (int iy = cy - ring; iy <= cy + ring; iy++) {
                if (iy < 0 || iy >= ny) continue;
                int onyface = (iy == cy - ring || iy == cy + ring);
                int zstep = (ix == xlo || ix == xhi || onyface) ? 1 : 2 * ring;
                if (zstep == 0) zstep = 1;
                for (int iz = cz - ring; iz <= cz + ring; iz += zstep) {
                    if (iz < 0 || iz >= nz) continue;
                    int c = ix * ny * nz + iy * nz + iz;
                    /* Squared-distance reject before the sqrt: an atom can
                       only improve `best` if its centre is within best+maxr,
                       and most atoms in a scanned cell are further than that.
                       The sqrt is the expensive part of this loop and it is
                       the whole search's hot path. */
                    double cut = best + maxr;
                    double cut2 = cut > 0 ? cut * cut : 0;
                    const double *ap = cellxyzr + 4 * cellstart[c];
                    for (int k = cellstart[c]; k < cellstart[c + 1]; k++, ap += 4) {
                        double dx = px - ap[0], dy = py - ap[1], dz = pz - ap[2];
                        double q = dx * dx + dy * dy + dz * dz;
                        if (cut > 0 && q > cut2) continue;
                        double d = sqrt(q) - ap[3];
                        if (d < best) { best = d; cut = best + maxr; cut2 = cut > 0 ? cut * cut : 0; }
                    }
                }
            }
        }
    }
    return best;
}

/* ---- in-plane Nelder-Mead maximizing clearance ------------------------ */
static double U[3], V[3], W[3];   /* orthonormal frame, W = channel vector */
static double O[3];               /* cpoint */

/* CAPSULE mode (HCAPEN, hcapen.f): the objective is not a point's clearance
 * but the largest-radius tube ("capsule") that fits along the SEGMENT from
 * the previous accepted centre (g_prevx/y/z) to the candidate point - i.e.
 * min over atoms of (perpendicular/endpoint distance to the segment, minus
 * vdw radius). raw = the physical capsule radius (what HOLE writes as the
 * .sph occupancy column); the search itself maximizes an "effective radius"
 * derived from it (equivalent-area circle of the swept stadium shape) only
 * when raw > 0 - a negative raw (segment already clipped by an atom) is
 * returned unchanged, exactly mirroring hcapen.f's own branch. */
static int g_capsule = 0;
static double g_prevx, g_prevy, g_prevz;

static double seg_clearance(double px, double py, double pz,
                            double bx, double by, double bz, double *out_dcent) {
    double ux = bx - px, uy = by - py, uz = bz - pz;
    double dcent = sqrt(ux*ux + uy*uy + uz*uz);
    *out_dcent = dcent;
    if (dcent < 1e-9) return clearance(px, py, pz);
    ux /= dcent; uy /= dcent; uz /= dcent;
    double mx = (px + bx) * 0.5, my = (py + by) * 0.5, mz = (pz + bz) * 0.5;
    int cx = (int)((mx - gx0) / cell), cy = (int)((my - gy0) / cell),
        cz = (int)((mz - gz0) / cell);
    double best = 1e30;
    double pad = 0.5 * dcent;      /* segment reaches this far past the midpoint */
    int mxx = (cx < 0 ? nx - cx : (cx >= nx ? cx : (cx > nx-1-cx ? cx : nx-1-cx)));
    int myy = (cy < 0 ? ny - cy : (cy >= ny ? cy : (cy > ny-1-cy ? cy : ny-1-cy)));
    int mzz = (cz < 0 ? nz - cz : (cz >= nz ? cz : (cz > nz-1-cz ? cz : nz-1-cz)));
    int maxring = (mxx > myy ? (mxx > mzz ? mxx : mzz) : (myy > mzz ? myy : mzz))
                + (int)(pad / cell) + 1;
    for (int ring = 0; ring <= maxring; ring++) {
        double lb = (ring - 1) * cell - maxr - pad;
        if (ring > 0 && lb > best) break;
        for (int ix = cx - ring; ix <= cx + ring; ix++) {
            if (ix < 0 || ix >= nx) continue;
            for (int iy = cy - ring; iy <= cy + ring; iy++) {
                if (iy < 0 || iy >= ny) continue;
                int onface = (ix == cx-ring || ix == cx+ring || iy == cy-ring || iy == cy+ring);
                int zstep = onface ? 1 : (ring > 0 ? 2*ring : 1);
                for (int iz = cz - ring; iz <= cz + ring; iz += zstep) {
                    if (iz < 0 || iz >= nz) continue;
                    int c = ix * ny * nz + iy * nz + iz;
                    for (int k = cellstart[c]; k < cellstart[c+1]; k++) {
                        int i = cellatom[k];
                        double rvx = ax[i]-px, rvy = ay[i]-py, rvz = az[i]-pz;
                        double rdotu = rvx*ux + rvy*uy + rvz*uz;
                        double d2;
                        if (rdotu < 0) {
                            d2 = rvx*rvx + rvy*rvy + rvz*rvz;
                        } else if (rdotu > dcent) {
                            double ex = bx-ax[i], ey = by-ay[i], ez = bz-az[i];
                            d2 = ex*ex + ey*ey + ez*ez;
                        } else {
                            double px2 = rvx - ux*rdotu, py2 = rvy - uy*rdotu, pz2 = rvz - uz*rdotu;
                            d2 = px2*px2 + py2*py2 + pz2*pz2;
                        }
                        double d = sqrt(d2) - ar[i];
                        if (d < best) best = d;
                    }
                }
            }
        }
    }
    return best;
}

static double f_at(double u, double v, double t) {
    double px = O[0] + u * U[0] + v * V[0] + t * W[0];
    double py = O[1] + u * U[1] + v * V[1] + t * W[1];
    double pz = O[2] + u * U[2] + v * V[2] + t * W[2];
    if (!g_capsule) return clearance(px, py, pz);
    double dcent;
    double raw = seg_clearance(px, py, pz, g_prevx, g_prevy, g_prevz, &dcent);
    if (raw <= 0.0) return raw;
    double area = M_PI * raw * raw + 2.0 * raw * dcent;
    return sqrt(area / M_PI);
}

static _Atomic long nm_evals = 0;
static double nm_cap = 1e30;   /* stop optimizing once clearance exceeds this
                                  (endrad + margin); outside the protein the
                                  objective is unbounded and expansion would
                                  otherwise run away */

/* maximize f in (u,v) at fixed t; returns best value, updates (u,v).
 * A two-pass "converge, then restart with a fresh simplex" variant was
 * tried and measured to never win in production (0/803 calls used its
 * second pass on the Nav bench, 0/594 on KcsA - the coarse-track ridge
 * probes and the sequential repair sweep already handle what a restart
 * was meant to catch). Removed rather than left dead: keeping unused
 * fallback logic around after superseding it is how "which algorithm is
 * actually running" stops being answerable from the code. */
/* How far a single call may move a vertex from where it started. nm_cap
   (checked below) stops the climb once its VALUE is clearly outside the
   pore, but a vestibule opening into bulk solvent can be near-flat in
   value over a wide area - every reflect/expand keeps finding a genuine,
   if marginal, improvement, so the value check never fires while the
   position runs away (measured on the Nav bench: 5 of 10 frames, up to a
   39A single-slice teleport, always the slice where the march happened to
   cross endrad). This is the classic Nelder-Mead pathology on a flat
   landscape, not a wall crossed - nm_connected() would pass. Every caller
   already treats this as a LOCAL refinement seeded at (*pu,*pv) (the
   coarse/fine track, the repair sweep, the offset probes reach further
   only by starting their OWN nm2d_once from an already-offset seed), so a
   hard cap on distance from that seed is the right invariant to enforce,
   not a looser one. */
#define NM_MAXSTEP 6.0
static double nm_maxstep2(void) {
    static double v = -1;
    if (v < 0) { const char *e = getenv("NM_MAXSTEP"); v = (e && atof(e) > 0) ? atof(e) : NM_MAXSTEP; v *= v; }
    return v;
}

static double nm2d_once(double *pu, double *pv, double t, double s) {
    double u0 = *pu, v0 = *pv;
    double px[3] = {*pu, *pu + s, *pu}, py[3] = {*pv, *pv, *pv + s}, pf[3];
    for (int i = 0; i < 3; i++) { pf[i] = f_at(px[i], py[i], t); nm_evals++; }
    double maxd2 = nm_maxstep2();
    for (int it = 0; it < 200; it++) {
        int hi = 0, lo = 0;
        for (int i = 1; i < 3; i++) {
            if (pf[i] < pf[hi]) hi = i;    /* worst = smallest clearance */
            if (pf[i] > pf[lo]) lo = i;    /* best */
        }
        if (pf[lo] > nm_cap) break;        /* out of the pore; caller stops */
        /* MEASURED, do not loosen without re-validating the FULL 10-frame
           corpus (not just the single frame + 2 static structures this
           was first checked against - that subset missed the failure).
           1e-3 looked free (identical accuracy on frame19/KcsA/9HNR, 16%
           faster) until it was run across all 10 real trajectory frames:
           frame 0 alone went from rms 0.027 to 1.34A. Root cause found,
           not just observed: the coarse pass's ridge-probe cascade tests
           candidates at fixed 1.5/3A offsets from the current seed, so a
           sub-0.01A difference in how precisely one early coarse point
           converged was enough to flip which quantised offset the probe
           cascade accepted - and that one flip put the rest of the march
           on a wholly different (still locally valid, still ascending)
           track through the vestibule. The probe mechanism amplifies
           convergence noise the same way CAPSULE mode's chained segment
           objective does (doc/07) - a small per-step perturbation, once
           accepted, compounds forward. 1e-4 does not trigger it on any
           of the 4 validated structures; treat any tolerance change as
           needing the same full-corpus check, not a 3-structure sample. */
        if (pf[lo] - pf[hi] < 1e-4) break;
        double cx = 0, cy = 0;
        for (int i = 0; i < 3; i++) if (i != hi) { cx += px[i]; cy += py[i]; }
        cx /= 2; cy /= 2;
        double rx = cx + (cx - px[hi]), ry = cy + (cy - py[hi]);
        double fr = f_at(rx, ry, t); nm_evals++;
        if ((rx-u0)*(rx-u0) + (ry-v0)*(ry-v0) > maxd2) fr = -1e30;
        if (fr > pf[lo]) {                                   /* expand */
            double ex = cx + 2 * (cx - px[hi]), ey = cy + 2 * (cy - py[hi]);
            double fe = f_at(ex, ey, t); nm_evals++;
            if ((ex-u0)*(ex-u0) + (ey-v0)*(ey-v0) > maxd2) fe = -1e30;
            if (fe > fr) { px[hi] = ex; py[hi] = ey; pf[hi] = fe; }
            else         { px[hi] = rx; py[hi] = ry; pf[hi] = fr; }
        } else if (fr > pf[hi]) {                            /* accept */
            px[hi] = rx; py[hi] = ry; pf[hi] = fr;
        } else {                                             /* contract */
            double kx = cx + 0.5 * (px[hi] - cx), ky = cy + 0.5 * (py[hi] - cy);
            double fk = f_at(kx, ky, t); nm_evals++;
            if (fk > pf[hi]) { px[hi] = kx; py[hi] = ky; pf[hi] = fk; }
            else {                                           /* shrink */
                for (int i = 0; i < 3; i++) {
                    if (i == lo) continue;
                    px[i] = px[lo] + 0.5 * (px[i] - px[lo]);
                    py[i] = py[lo] + 0.5 * (py[i] - py[lo]);
                    pf[i] = f_at(px[i], py[i], t); nm_evals++;
                }
            }
        }
    }
    int lo = 0;
    for (int i = 1; i < 3; i++) if (pf[i] > pf[lo]) lo = i;
    *pu = px[lo]; *pv = py[lo];
    return pf[lo];
}

static double now_ms(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1e3 + ts.tv_nsec / 1e6;
}

/* ---- reusable pipeline pieces ---------------------------------------- */

/* orthonormal frame from cpoint + cvect, then the atom grid */
static void nm_setup_frame(const double cp[3], const double cv[3]) {
    O[0] = cp[0]; O[1] = cp[1]; O[2] = cp[2];
    W[0] = cv[0]; W[1] = cv[1]; W[2] = cv[2];
    double wn = sqrt(W[0]*W[0] + W[1]*W[1] + W[2]*W[2]);
    for (int i = 0; i < 3; i++) W[i] /= wn;
    double t[3] = {1, 0, 0};                /* any vector not parallel to W */
    if (fabs(W[0]) > 0.9) { t[0] = 0; t[1] = 1; }
    U[0] = W[1]*t[2] - W[2]*t[1]; U[1] = W[2]*t[0] - W[0]*t[2]; U[2] = W[0]*t[1] - W[1]*t[0];
    double un = sqrt(U[0]*U[0] + U[1]*U[1] + U[2]*U[2]);
    for (int i = 0; i < 3; i++) U[i] /= un;
    V[0] = W[1]*U[2] - W[2]*U[1]; V[1] = W[2]*U[0] - W[0]*U[2]; V[2] = W[0]*U[1] - W[1]*U[0];
}

/* march results */
#define MAXSL 8192
static double slt[MAXSL], slx[MAXSL], sly[MAXSL], slz[MAXSL], slr[MAXSL];
static int nslice = 0;

/* Two-phase march (fills sl*): clipstop >= endrad lets the march continue
 * past the endrad crossing so the caller gets HOLE-style oversized "clip"
 * spheres for the mouths; clipstop == endrad reproduces the plain profile. */
static double tm_coarse, tm_fine, tm_repair;   /* stage timings, ms */

/* Same passage or a different one? Two centres in one slice plane belong to
   the same pore if the straight segment between them never enters an atom:
   clearance stays positive the whole way. A better clearance across a WALL is
   a neighbouring passage, and following it is how a track abandons the pore
   it was asked to trace (frame 90 of the Nav trajectory: a 0.23 A constriction
   replaced by a 1.12 A one in a different channel). Sampled every 0.25 A. */
static int nm_connected(double u1, double v1, double u2, double v2, double t) {
    double du = u2 - u1, dv = v2 - v1, len = sqrt(du*du + dv*dv);
    int n = (int)(len / 0.25) + 1;
    for (int i = 0; i <= n; i++) {
        double s = (double)i / n;
        if (f_at(u1 + s*du, v1 + s*dv, t) <= 0.0) return 0;
    }
    return 1;
}

static void nm_march(double sample, double endrad, double clipstop) {
    double _t0 = now_ms();
    nm_cap = clipstop + 2.0;
    nslice = 0;
    /* Phase A (sequential, cheap): coarse march both directions from cpoint,
     * seeded from the previous coarse slice, to lay down approximate centres.
     * Phase B (parallel): every fine slice refined independently, seeded by
     * interpolating the coarse centres at its t. Slices are independent once
     * seeded, so this is an embarrassingly parallel loop. */
    static double ct[MAXSL], cu[MAXSL], cv[MAXSL];   /* coarse centres */
    int ncoarse = 0;
    double csamp = sample * 8.0;
    { const char *e = getenv("NM_CSAMP_MULT"); if (e && atof(e) > 0) csamp = sample * atof(e); }
    if (csamp < 1.0) csamp = 1.0;
    /* NOTE (measured, do not re-add): seeding frames 2+ from the previous
       frame's track instead of a fresh coarse pass was tried and is a NET
       LOSS - inherited ridge errors accumulate (f06 rms 0.18 -> 0.59), and
       any end-margin compounds frame over frame (z span grew 40 -> 52 A in
       10 frames). Per-frame independence wins on accuracy AND simplicity. */
    for (int dir = -1; dir <= 1; dir += 2) {
        double u = 0, v = 0;
        for (int k = (dir == 1 ? 0 : 1); k < 2000; k++) {
            double tt = dir * k * csamp;
            double r = nm2d_once(&u, &v, tt, 0.4);
            /* Local multi-start, coarse pass only: HOLE's annealing crosses
               small ridges and lands on slightly wider centres a couple of A
               away (measured on the Nav trajectory: r deficits of 0.1-0.3 A
               track-wide when NM stays on the lesser ridge). Offsets stay at
               pore-interior scale — large offsets jump to bulk at the mouths
               (see the note in the fine pass). The winning track propagates
               to every fine slice through the seeds. */
            /* Offset scale and round count are tunable so the cascade can be
               swept against HOLE rather than assumed; defaults are the values
               validated against HOLE on the benchmark structures. */
            double _sc = 1.0; int _rounds = 2;
            { const char *e = getenv("NM_PROBE_SCALE"); if (e && atof(e) > 0) _sc = atof(e);
              const char *g = getenv("NM_PROBE_ROUNDS"); if (g && atoi(g) > 0) _rounds = atoi(g); }
            const double mo[8][2] = {{1.5*_sc,0},{-1.5*_sc,0},{0,1.5*_sc},{0,-1.5*_sc},
                                     {3*_sc,0},{-3*_sc,0},{0,3*_sc},{0,-3*_sc}};
            if (r <= 0.55 * endrad) {      /* core only */
                /* two parallel probe rounds, re-centred on each round's
                   winner: the cascade is what walks the track across a
                   clearance ridge, and running each round's 8 single-pass
                   probes on the thread pool makes a round cost ~1 probe */
                for (int round = 0; round < _rounds; round++) {
                    double bu[8], bv[8], br[8];
                    #pragma omp parallel for schedule(static)
                    for (int q = 0; q < 8; q++) {
                        bu[q] = u + mo[q][0]; bv[q] = v + mo[q][1];
                        /* one clearance eval rejects offsets starting inside
                           the wall - a simplex started there burns its whole
                           budget crawling out, and a real centre never sits
                           at negative clearance */
                        if (f_at(bu[q], bv[q], tt) < 0.3) { br[q] = -1e30; continue; }
                        br[q] = nm2d_once(&bu[q], &bv[q], tt, 0.4);
                    }
                    int best = -1;
                    for (int q = 0; q < 8; q++)
                        if (br[q] > r && br[q] <= nm_cap) { r = br[q]; best = q; }
                    if (best < 0) break;
                    u = bu[best]; v = bv[best];
                }
            }
            ct[ncoarse] = tt; cu[ncoarse] = u; cv[ncoarse] = v; ncoarse++;
            if (getenv("NM_DEBUG_COARSE"))
                fprintf(stderr, "coarse t=%7.2f u=%7.2f v=%7.2f r=%.3f\n", tt, u, v, r);
            if (r > clipstop || ncoarse >= MAXSL) break;
        }
    }
    tm_coarse = now_ms() - _t0; _t0 = now_ms();
    /* sort coarse centres by t (two small runs; insertion sort is fine) */
    for (int i = 1; i < ncoarse; i++) {
        double a = ct[i], b = cu[i], c = cv[i]; int j = i - 1;
        while (j >= 0 && ct[j] > a) { ct[j+1]=ct[j]; cu[j+1]=cu[j]; cv[j+1]=cv[j]; j--; }
        ct[j+1] = a; cu[j+1] = b; cv[j+1] = c;
    }
    double tmin = ct[0], tmax = ct[ncoarse - 1];
    int k0 = (int)floor(tmin / sample), k1 = (int)ceil(tmax / sample);
    int nfine = k1 - k0 + 1;
    static double fu[MAXSL], fv[MAXSL], fr[MAXSL];
    if (nfine > MAXSL) nfine = MAXSL;
    if (getenv("NM_SEQ") && atoi(getenv("NM_SEQ"))) {
        /* Sequential fine pass: each slice seeded from the previous FINE
           slice, the coarse pass's probe cascade applied at every narrow
           slice, the two march directions on two threads. This is the
           algorithm's own definition (doc/03) and the analogue of HOLE's
           limited-mobility annealing: the track can only move gradually,
           so a basin that swaps along t is followed (frame 60) and a wider
           pocket contiguous in one plane is not jumped into (frame 90).
           No repair sweep: continuity holds by construction. */
        int kc = -k0;   /* index of t = 0 */
        #pragma omp parallel for num_threads(2) schedule(static)
        for (int dir = 0; dir < 2; dir++) {
            int step = dir ? 1 : -1;
            double u = 0, v = 0;
            for (int k = kc; k >= 0 && k < nfine; k += step) {
                double tt = (k0 + k) * sample;
                double r = nm2d_once(&u, &v, tt, 0.4);
                int _sp = 1; { const char *e = getenv("NM_SEQ_PROBE"); if (e) _sp = atoi(e); }
                if (_sp && r <= 0.55 * endrad) {
                    static const double mo[8][2] = {{1.5,0},{-1.5,0},{0,1.5},{0,-1.5},
                                                    {3,0},{-3,0},{0,3},{0,-3}};
                    for (int round = 0; round < 2; round++) {
                        double bu = u, bv = v, br = r;
                        for (int q = 0; q < 8; q++) {
                            double pu = u + mo[q][0], pv = v + mo[q][1];
                            if (f_at(pu, pv, tt) < 0.3) continue;
                            double pr = nm2d_once(&pu, &pv, tt, 0.4);
                            if (pr > br) { br = pr; bu = pu; bv = pv; }
                        }
                        if (br <= r) break;
                        u = bu; v = bv; r = br;
                    }
                }
                fr[k] = r; fu[k] = u; fv[k] = v;
                if (dir == 0 && k == kc) { /* t=0 shared: the forward march re-seeds from it */ }
            }
        }
        tm_fine = now_ms() - _t0; _t0 = now_ms();
        nslice = 0;
        for (int k = 0; k < nfine; k++) {
            if (nslice >= MAXSL) break;
            double tt = (k0 + k) * sample;
            slt[nslice] = tt;
            slx[nslice] = O[0] + fu[k]*U[0] + fv[k]*V[0] + tt*W[0];
            sly[nslice] = O[1] + fu[k]*U[1] + fv[k]*V[1] + tt*W[1];
            slz[nslice] = O[2] + fu[k]*U[2] + fv[k]*V[2] + tt*W[2];
            slr[nslice] = fr[k];
            nslice++;
        }
        tm_repair = 0;
        return;
    }
    #pragma omp parallel for schedule(dynamic, 8)
    for (int k = 0; k < nfine; k++) {
        double tt = (k0 + k) * sample;
        /* seed from coarse centres: nearest-t linear interpolation */
        int j = 0;
        while (j < ncoarse - 2 && ct[j + 1] < tt) j++;
        double span = ct[j + 1] - ct[j];
        double w2 = span > 0 ? (tt - ct[j]) / span : 0;
        if (w2 < 0) w2 = 0; if (w2 > 1) w2 = 1;
        double u = cu[j] + w2 * (cu[j + 1] - cu[j]);
        double v = cv[j] + w2 * (cv[j + 1] - cv[j]);
        fr[k] = nm2d_once(&u, &v, tt, 0.4);
        fu[k] = u; fv[k] = v;
        /* Constriction-slice probe, at the coarse cascade's own offsets.
           The coarse pass samples t every csamp and its probes pick the best
           basin THERE; when the better basin only exists between two coarse
           samples every fine slice in that stretch inherits the lesser one and
           a single seeded simplex cannot leave it. Measured on frame 60 of the
           Nav trajectory: the port sat at 0.767 A while a basin of 1.43 A lay
           3.25 A away across a 0.51 A saddle - no wall, just out of reach of a
           local step. A 0.75 A probe here was tried first and made things
           worse; the offsets must be the ones that reach the next basin.
           Gated to narrow slices, where a 3 A step cannot reach bulk (the
           note below is about wide slices, where it can). */
        {
            double _pg = 0.0, _ps = 1.0;   /* off: see VALIDATION.md, it hops into side pockets */
            { const char *e = getenv("NM_FINE_GATE");  if (e && atof(e) >= 0) _pg = atof(e);
              const char *g = getenv("NM_FINE_SCALE"); if (g && atof(g) >  0) _ps = atof(g); }
            if (_pg > 0 && fr[k] < _pg * endrad) {
                const double fo[8][2] = {{1.5*_ps,0},{-1.5*_ps,0},{0,1.5*_ps},{0,-1.5*_ps},
                                         {3*_ps,0},{-3*_ps,0},{0,3*_ps},{0,-3*_ps}};
                for (int q = 0; q < 8; q++) {
                    double pu = u + fo[q][0], pv = v + fo[q][1];
                    if (f_at(pu, pv, tt) < 0.3) continue;
                    double pr = nm2d_once(&pu, &pv, tt, 0.4);
                    if (pr > fr[k] && nm_connected(u, v, pu, pv, tt)) { fr[k] = pr; fu[k] = pu; fv[k] = pv; }
                }
            }
        }
        /* NOTE (measured, do not "fix"): no multi-start here. At vestibule
           slices the widest CONNECTED sphere in the plane sits in bulk -
           the mouth is open to bulk by definition - so any global search
           hijacks the slice (r jumps 8.7 -> 19.4 at z=34 on the Nav bench,
           truncating the march). HOLE only avoids this because annealing is
           a limited-mobility local walk; seeded NM is the same by
           construction, which is the behaviour we want. */
    }
    tm_fine = now_ms() - _t0; _t0 = now_ms();
    /* Sequential repair sweep. The algorithm's definition is "seeded from
       the previous slice"; the parallel pass approximates that with coarse
       interpolated seeds, which occasionally lands a slice in a different
       basin (measured: track jumps up to 35 A for a few slices on 3 of 10
       trajectory frames). Walk outward from t=0 and re-run any slice whose
       centre jumps too far from its inward neighbour, seeded from that
       neighbour - restoring sequential semantics exactly where they
       diverged, at ~zero cost on frames that already agree. */
    {
        int kc0 = -k0;
        double jumptol = 1.0;
        { const char *e = getenv("NM_JUMPTOL"); if (e && atof(e) > 0) jumptol = atof(e); }
        for (int dir2 = -1; dir2 <= 1; dir2 += 2) {
            int prev = kc0;
            for (int k = kc0 + dir2; k >= 0 && k < nfine; k += dir2) {
                double du = fu[k] - fu[prev], dv = fv[k] - fv[prev];
                if (du*du + dv*dv > jumptol*jumptol) {
                    double u = fu[prev], v = fv[prev];
                    double rr = nm2d_once(&u, &v, (k0 + k) * sample, 0.4);
                    /* At a NARROW slice keep whichever of the two is the
                       better clearance, the same criterion everything else
                       uses. Replacing unconditionally let one poor re-run
                       poison every slice outward of it: each next slice then
                       "jumped" relative to its now-poor neighbour, was re-run
                       from it, and inherited the same basin (frame 60: a
                       stretch that had 1.03-1.24 A was dragged to 0.77 A).
                       A WIDE slice keeps the replacement: there the jump the
                       sweep exists for is a hijack into bulk, whose clearance
                       is larger by definition, so "better" would keep it. */
                    /* Unconditional replace, on purpose. "Keep the slice's own
                       result when its clearance is higher" was tried, with and
                       without an in-plane connectivity test: it lifts frame 60
                       (0.77 -> 0.99 A) but hops into contiguous side pockets on
                       frames 70 and 90 (0.61 -> 0.98, 0.23 -> 0.61 A). Continuity
                       along t is the property that defines the track. */
                    fr[k] = rr; fu[k] = u; fv[k] = v;
                }
                prev = k;
            }
        }
    }
    tm_repair = now_ms() - _t0;
    /* trim to the contiguous pore: walk outward from t=0, past endrad keep
       going while r <= clipstop (those become mouth clip spheres), stop at
       the first slice beyond clipstop (keep the crossing slice) */
    int kc = -k0;                      /* index of t = 0 */
    int lo = kc, hi = kc;
    while (lo > 0 && fr[lo] <= clipstop) lo--;
    while (hi < nfine - 1 && fr[hi] <= clipstop) hi++;
    for (int k = lo; k <= hi && nslice < MAXSL; k++) {
        double tt = (k0 + k) * sample;
        slt[nslice] = tt;
        slx[nslice] = O[0] + fu[k]*U[0] + fv[k]*V[0] + tt*W[0];
        sly[nslice] = O[1] + fu[k]*U[1] + fv[k]*V[1] + tt*W[1];
        slz[nslice] = O[2] + fu[k]*U[2] + fv[k]*V[2] + tt*W[2];
        slr[nslice] = fr[k];
        nslice++;
    }
}

#ifndef NEWHOLE_EMBED
#include "nm_holeout.h"

static void usage(const char *a0) {
    fprintf(stderr,
        "usage: %s pdb rad cx cy cz vx vy vz sample endrad [options]\n"
        "  --sph FILE        write HOLE-format hole_out.sph\n"
        "  --tsv FILE        write the 8-column hole_profile.tsv\n"
        "  --conn PROBE GRID Connolly pass (GRID 0 = 0.7*PROBE), needs --sph/--tsv\n"
        "  --ignore R1,R2    drop residues by name before the search (HOLE IGNORE)\n"
        "  --centres FILE    skip the search, take 't x y z r' centres from FILE\n"
        "  --quiet           no slice listing on stdout\n", a0);
}

#ifdef VMDPATHFINDER_MULTICALL
int nm_search_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
    if (argc < 11) { usage(argv[0]); return 2; }
    nm_set_wait_policy(argv);
    const char *sph = NULL, *tsv = NULL, *ignore = NULL, *centres = NULL;
    double probe = 0, grid = 0; int conn = 0, quiet = 0;
    for (int a = 11; a < argc; a++) {
        if (!strcmp(argv[a], "--sph") && a + 1 < argc) sph = argv[++a];
        else if (!strcmp(argv[a], "--tsv") && a + 1 < argc) tsv = argv[++a];
        else if (!strcmp(argv[a], "--ignore") && a + 1 < argc) ignore = argv[++a];
        else if (!strcmp(argv[a], "--centres") && a + 1 < argc) centres = argv[++a];
        else if (!strcmp(argv[a], "--conn") && a + 2 < argc) { conn = 1; probe = atof(argv[++a]); grid = atof(argv[++a]); }
        else if (!strcmp(argv[a], "--quiet")) quiet = 1;
        else { usage(argv[0]); return 2; }
    }
    if (conn) {
        if (probe <= 0) probe = 1.15;
        if (grid <= 0) grid = 0.7 * probe;
    }
    nm_set_default_threads();
    double t0 = now_ms();
    load_rad(argv[2]);
    load_pdb(argv[1]);
    int dropped = ignore ? ho_apply_ignore(ignore, aresn) : 0;
    double t1 = now_ms();
    double cp[3] = {atof(argv[3]), atof(argv[4]), atof(argv[5])};
    double cv[3] = {atof(argv[6]), atof(argv[7]), atof(argv[8])};
    double sample = atof(argv[9]), endrad = atof(argv[10]);
    nm_setup_frame(cp, cv);
    build_grid();
    double t2 = now_ms();
    if (centres) { if (!ho_load_centres(centres)) return 1; }
    else nm_march(sample, endrad, endrad);
    double t3 = now_ms();
    if (nslice < 1) { fprintf(stderr, "no slices found\n"); return 1; }
    ho_layout L; ho_layout_init(&L, sample);
    if (conn) ho_run_conn(&L, endrad, probe, grid);
    double t4 = now_ms();
    if (sph && !ho_write_sph(sph, &L, endrad)) return 1;
    if (tsv && !ho_write_tsv(tsv, &L, endrad)) return 1;
    if (!quiet)
        for (int i = 0; i < nslice; i++)
            printf("%.3f %.5f %.5f %.5f %.5f\n", slt[i], slx[i], sly[i], slz[i], slr[i]);
    fprintf(stderr, "atoms %d (ignored %d)  load+radii %.1f ms  grid %.1f ms  search %.1f ms  (%ld evals)",
            natom, dropped, t1 - t0, t2 - t1, t3 - t2, (long)nm_evals);
    if (conn) fprintf(stderr, "  conn %.1f ms (%d slices)", t4 - t3, nslice);
    fprintf(stderr, "\n");
    ho_layout_free(&L);
    return 0;
}
#endif
