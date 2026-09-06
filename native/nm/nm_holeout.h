/* nm_holeout.h - HOLE-format output for nm_search, plus a literal port of
 * HOLE's per-slice Connolly step.
 *
 * Included by nm_search.c after clearance()/the slice arrays exist. Three
 * things live here:
 *
 *   1. write_sph()  - HOLE's hole_out.sph layout (holcal.f order): +ve slices,
 *      ADDEND's escape spheres (resid -888), the start slice again, -ve
 *      slices, ADDEND again. Every consumer in the plugin (sph_process, the
 *      lobe classifier, Ion Flow) reads this layout, so it is reproduced
 *      exactly rather than "close enough".
 *   2. write_tsv()  - the 8-column hole_profile.tsv the plugin parses when no
 *      hole_out.txt exists (parse_profile_from_tsv): slice rows plus HOLE's
 *      mid-point rows, cen_line_d arc length, sum{s/area}, and the three
 *      Connolly columns when a conn pass ran.
 *   3. concal()/coarea()/addend() - ported line for line from the plugin's
 *      bit-exact Tcl port of concal.f / coarea.f / addend.f, including
 *      dist2_buggy() (HOLE adds dz+dz instead of dz*dz in one duplicate test;
 *      keeping the bug is what makes the dot sets comparable). The Connolly
 *      pass is the only part that is parallel: slices are independent, so
 *      OpenMP runs them concurrently and the writer serialises afterwards.
 *
 * Conventions shared with the rest of nm_search: slt[] is the axial offset
 * from CPOINT along the unit channel vector W, slx/sly/slz the sphere centre,
 * slr its radius. HOLE's "coord" column is the absolute projection
 * centre.W, not slt, so the TSV uses that.
 */

#include <ctype.h>

/* ------------------------------------------------------------------ frame */

/* calper.f: HOLE's own perpendicular pair for a unit channel vector. The
   Connolly grid walks along these two vectors, so the port must use HOLE's
   choice of basis, not nm_search's U/V. */
static void ho_calper(const double v[3], double e[3], double n[3]) {
    if (fabs(v[1]) < 1e-6 && fabs(v[2]) < 1e-6) { e[0] = 0; e[1] = 1; e[2] = 0; }
    else                                        { e[0] = 1; e[1] = 0; e[2] = 0; }
    double cdotp = v[0]*e[0] + v[1]*e[1] + v[2]*e[2];
    for (int i = 0; i < 3; i++) e[i] -= cdotp * v[i];
    double en = sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    for (int i = 0; i < 3; i++) e[i] /= en;
    n[0] = e[1]*v[2] - e[2]*v[1];
    n[1] = e[2]*v[0] - e[0]*v[2];
    n[2] = e[0]*v[1] - e[1]*v[0];
}

static double ho_dist2_buggy(double ax_, double ay_, double az_,
                             double bx_, double by_, double bz_) {
    double d1 = ax_ - bx_, d2 = ay_ - by_, d3 = az_ - bz_;
    return d1*d1 + d2*d2 + d3 + d3;
}

static void ho_sph_line(FILE *f, int irec, double x, double y, double z,
                        double rad1, double rad2) {
    fprintf(f, "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
            1, "QSS", "SPH", "S", irec, x, y, z, rad1, rad2);
}

/* ---------------------------------------------------------------- concal */

#define HO_SCOMAX 10000

/* Duplicate test, exact but not O(n): the buggy metric is d1*d1 + d2*d2 + 2*d3
   with d3 = z_new - z_s, so a point s can only pass the "< thr" test when its
   world-xy distance satisfies d1*d1 + d2*d2 < thr - 2*(z_new - zmax). Points
   outside that xy radius are skipped; the ones inside get the original test. */
typedef struct { int *head, *next; double h, zmax; } ho_hash;
#define HO_HBITS 16
static unsigned ho_hkey(int ix, int iy) {
    return ((unsigned)ix * 73856093u ^ (unsigned)iy * 19349663u) & ((1u << HO_HBITS) - 1);
}
static void ho_hash_init(ho_hash *H, double h) {
    H->head = malloc((1 << HO_HBITS) * sizeof(int));
    for (int i = 0; i < (1 << HO_HBITS); i++) H->head[i] = -1;
    H->next = malloc(HO_SCOMAX * sizeof(int));
    H->h = h; H->zmax = -1e300;
}
static void ho_hash_free(ho_hash *H) { free(H->head); free(H->next); }
static void ho_hash_add(ho_hash *H, int idx, double x, double y, double z) {
    unsigned k = ho_hkey((int)floor(x / H->h), (int)floor(y / H->h));
    H->next[idx] = H->head[k]; H->head[k] = idx;
    if (z > H->zmax) H->zmax = z;
}
static int ho_dup(const ho_hash *H, const double *px, const double *py, const double *pz,
                  double x, double y, double z, double thr) {
    double r2 = thr - 2.0 * (z - H->zmax);
    if (r2 <= 0) return 0;
    double r = sqrt(r2) + 1e-9;
    int ix0 = (int)floor((x - r) / H->h), ix1 = (int)floor((x + r) / H->h);
    int iy0 = (int)floor((y - r) / H->h), iy1 = (int)floor((y + r) / H->h);
    for (int ix = ix0; ix <= ix1; ix++)
        for (int iy = iy0; iy <= iy1; iy++)
            for (int s = H->head[ho_hkey(ix, iy)]; s >= 0; s = H->next[s])
                if (ho_dist2_buggy(x, y, z, px[s], py[s], pz[s]) < thr) return 1;
    return 0;
}

typedef struct {
    int n;                 /* points found (0 when fallback) */
    double *px, *py, *pz, *pr;
    unsigned char *act;    /* coarea's final active flags */
    double requiv;         /* 1e6 = escaped, 0 = no circles */
    double rad0;
    int fallback, escaped;
} ho_slice;

static void ho_slice_free(ho_slice *s) {
    free(s->px); free(s->py); free(s->pz); free(s->pr); free(s->act);
    memset(s, 0, sizeof *s);
}

/* concal.f: grow the Connolly point set outward from the sphere centre in the
   slice plane. Returns 0 on fallback (centre clearance below the probe). */
static int ho_concal(double cx, double cy, double cz, const double e[3],
                     const double nn[3], double endrad, double probe,
                     double grid, ho_slice *out) {
    memset(out, 0, sizeof *out);
    double endrp3 = endrad + 3.0;
    double rad0 = clearance(cx, cy, cz);
    out->rad0 = rad0;
    if (rad0 < probe) { out->fallback = 1; out->requiv = rad0; return 0; }
    double *px = malloc(HO_SCOMAX * sizeof(double)), *py = malloc(HO_SCOMAX * sizeof(double)),
           *pz = malloc(HO_SCOMAX * sizeof(double)), *pr = malloc(HO_SCOMAX * sizeof(double));
    unsigned char *pdo = malloc(HO_SCOMAX), *act = malloc(HO_SCOMAX);
    int sconum = 1, escaped = 0, sconxt = 0;
    ho_hash H; ho_hash_init(&H, grid > 0.5 ? grid : 0.5);
    px[0] = cx; py[0] = cy; pz[0] = cz; pr[0] = rad0; pdo[0] = 1;
    ho_hash_add(&H, 0, cx, cy, cz);
    for (;;) {
        double r = pr[sconxt];
        if (r < endrp3 && pdo[sconxt]) {
            pdo[sconxt] = 0;
            double bx = px[sconxt], by = py[sconxt], bz = pz[sconxt];
            for (int ncount = -1; ncount <= 1 && !escaped; ncount++) {
                for (int ecount = -1; ecount <= 1 && !escaped; ecount++) {
                    double tx = bx + ecount*grid*e[0] + ncount*grid*nn[0];
                    double ty = by + ecount*grid*e[1] + ncount*grid*nn[1];
                    double tz = bz + ecount*grid*e[2] + ncount*grid*nn[2];
                    if (ho_dup(&H, px, py, pz, tx, ty, tz, grid/1000.0)) continue;
                    double newrad = clearance(tx, ty, tz);
                    if (newrad >= probe) {
                        if (sconum + 1 > HO_SCOMAX) goto full;
                        px[sconum] = tx; py[sconum] = ty; pz[sconum] = tz;
                        pr[sconum] = newrad; pdo[sconum] = 1; ho_hash_add(&H, sconum, tx, ty, tz); sconum++;
                        if (newrad > endrp3) escaped = 1;
                    } else {
                        double vxu = tx - bx, vyu = ty - by, vzu = tz - bz;
                        double vn = sqrt(vxu*vxu + vyu*vyu + vzu*vzu);
                        vxu /= vn; vyu /= vn; vzu /= vn;
                        double delta = 0.25 * grid, excess = 0.0;
                        double sx = tx, sy = ty, sz = tz;
                        for (int cc = 0; cc < 100; cc++) {
                            sx = bx + delta*vxu; sy = by + delta*vyu; sz = bz + delta*vzu;
                            double srad = clearance(sx, sy, sz);
                            excess = srad - (probe + 0.0001);
                            if (fabs(excess) < 0.0005) break;
                            delta += excess;
                        }
                        if (fabs(excess) < 0.0005 && delta > 0.0) {
                            if (!ho_dup(&H, px, py, pz, sx, sy, sz, 0.09)) {
                                if (sconum + 1 > HO_SCOMAX) goto full;
                                px[sconum] = sx; py[sconum] = sy; pz[sconum] = sz;
                                pr[sconum] = excess + probe + 0.0001; pdo[sconum] = 0;
                                ho_hash_add(&H, sconum, sx, sy, sz); sconum++;
                            }
                        }
                    }
                }
            }
        }
        if (escaped) break;
        double maxrad = -1e10; int nxt = -1;
        for (int s = 0; s < sconum; s++)
            if (pdo[s] && pr[s] > maxrad) { maxrad = pr[s]; nxt = s; }
        if (nxt < 0) break;
        if (pr[nxt] >= endrp3) { escaped = 1; break; }   /* only the centre itself can be
                                                           this wide (any other such point
                                                           set escaped on insertion); HOLE's
                                                           search never hands concal one */
        sconxt = nxt;
    }
full:
    ho_hash_free(&H);
    for (int s = 0; s < sconum; s++) act[s] = 1;
    int last = sconum - 1;
    if (pr[last] > endrp3) {
        double delim2 = pr[last] * pr[last];
        for (int s = 0; s < last; s++)
            if (ho_dist2_buggy(px[last], py[last], pz[last], px[s], py[s], pz[s]) < delim2) act[s] = 0;
    }
    if (pr[0] > 2.5 * probe) {
        double d = pr[0] - 1.5 * probe, delim2 = d * d;
        for (int s = 1; s < sconum; s++) {
            if (!act[s]) continue;
            if (ho_dist2_buggy(px[0], py[0], pz[0], px[s], py[s], pz[s]) < delim2) act[s] = 0;
        }
    }
    free(pdo);
    out->n = sconum; out->px = px; out->py = py; out->pz = pz; out->pr = pr; out->act = act;
    out->escaped = escaped;
    return 1;
}

/* ---------------------------------------------------------------- coarea */

#define HO_SNMAX 90000

/* coarea.f: Requiv from the area of the union of the in-plane circles,
   integrated on a 0.25 A grid with adaptive refinement of the grey cells.
   Updates s->act (points whose circle leaves the plane are dropped) and
   sets s->requiv.

   Same arithmetic as the Fortran/Tcl original at every cell, reached faster:
   the first raster is a recursive block fill (a block wholly inside one
   shrunken circle is black, a block no widened circle reaches is white,
   anything else splits down to single cells that run the original per-cell
   loop over the circles that can reach them), and each grey cell keeps that
   candidate list for the refinement cycles. The block shortcuts are only
   taken when they imply the per-cell result with a margin, so the black/grey
   classification, the area sums (same count of identical terms) and the
   grey-list order (raster order) are the ones the plain loop produces. */

typedef struct {
    int nv; const int *vi;                 /* valid circle indices, ascending */
    const double *ec, *nc, *cr;
    const double *edc, *ndc;               /* cell origins, by repeated addition */
    double size1, m;                       /* m = 0.5*root2*size1 */
    unsigned char *mask;                   /* per cell: 0 white 1 grey 2 black */
    int tripN;
    int *cand; int *cand_off, *cand_n;     /* per-cell candidate lists (grey cells) */
    int ncand, capcand;
} ho_area;

static void ho_cand_push(ho_area *A, int cell, const int *sub, int nsub) {
    if (A->ncand + nsub > A->capcand) {
        A->capcand = 2 * (A->ncand + nsub) + 1024;
        A->cand = realloc(A->cand, A->capcand * sizeof(int));
    }
    A->cand_off[cell] = A->ncand; A->cand_n[cell] = nsub;
    memcpy(A->cand + A->ncand, sub, nsub * sizeof(int));
    A->ncand += nsub;
}

static void ho_block(ho_area *A, int ie0, int ie1, int jn0, int jn1, const int *cand, int ncand) {
    double s1 = A->size1, m = A->m;
    double e_lo = A->edc[ie0] + 0.5*s1, e_hi = A->edc[ie1-1] + 0.5*s1;
    double n_lo = A->ndc[jn0] + 0.5*s1, n_hi = A->ndc[jn1-1] + 0.5*s1;
    int *sub = malloc((ncand > 0 ? ncand : 1) * sizeof(int)); int nsub = 0, fullblack = 0;
    for (int q = 0; q < ncand; q++) {
        int i = cand[q];
        double ce = A->ec[i], cn = A->nc[i], cr = A->cr[i];
        double dme = ce < e_lo ? e_lo - ce : (ce > e_hi ? ce - e_hi : 0.0);
        double dmn = cn < n_lo ? n_lo - cn : (cn > n_hi ? cn - n_hi : 0.0);
        double xe = fabs(e_lo - ce) > fabs(e_hi - ce) ? fabs(e_lo - ce) : fabs(e_hi - ce);
        double xn = fabs(n_lo - cn) > fabs(n_hi - cn) ? fabs(n_lo - cn) : fabs(n_hi - cn);
        double dmin = sqrt(dme*dme + dmn*dmn), dmax = sqrt(xe*xe + xn*xn);
        /* reach = widened circle plus one cell diagonal, so the leaf list also
           covers every sub-cell centre of the refinement */
        if (dmin - 1e-7 >= cr + m + m) continue;
        sub[nsub++] = i;
        if (dmax + 1e-7 < fabs(cr - m)) { fullblack = 1; break; }
    }
    if (fullblack) {
        for (int ie = ie0; ie < ie1; ie++)
            for (int jn = jn0; jn < jn1; jn++) A->mask[ie * A->tripN + jn] = 2;
        free(sub);
        return;
    }
    if (nsub == 0) { free(sub); return; }
    if (ie1 - ie0 == 1 && jn1 - jn0 == 1) {
        double cenE = A->edc[ie0] + 0.5*s1, cenN = A->ndc[jn0] + 0.5*s1;
        int lblack = 0, isblack = 0;
        for (int q = 0; q < nsub; q++) {
            int i = sub[q];
            double cr = A->cr[i], de = cenE - A->ec[i], dn = cenN - A->nc[i];
            double dist2 = de*de + dn*dn, w = cr - m;
            if (dist2 < w*w) { isblack = 1; break; }
            w = cr + m;
            if (dist2 < w*w) lblack = 1;
        }
        int cell = ie0 * A->tripN + jn0;
        if (isblack) A->mask[cell] = 2;
        else if (lblack) { A->mask[cell] = 1; ho_cand_push(A, cell, sub, nsub); }
        free(sub);
        return;
    }
    if (ie1 - ie0 >= jn1 - jn0) {
        int mid = (ie0 + ie1) / 2;
        ho_block(A, ie0, mid, jn0, jn1, sub, nsub);
        ho_block(A, mid, ie1, jn0, jn1, sub, nsub);
    } else {
        int mid = (jn0 + jn1) / 2;
        ho_block(A, ie0, ie1, jn0, mid, sub, nsub);
        ho_block(A, ie0, ie1, mid, jn1, sub, nsub);
    }
    free(sub);
}

typedef struct { double bx, by; int root; } ho_cell;

static void ho_coarea(ho_slice *s, double cx, double cy, double cz,
                      const double e[3], const double nn[3], const double v[3],
                      double sample, double endrad) {
    int n = s->n;
    const double pi = 2.0 * asin(1.0), root2 = sqrt(2.0);
    double *circrad = calloc(n, sizeof(double)), *ec = calloc(n, sizeof(double)),
           *nc = calloc(n, sizeof(double));
    unsigned char *act = s->act;
    for (int i = 0; i < n; i++) {
        if (!act[i]) continue;
        double tx = s->px[i] - cx, ty = s->py[i] - cy, tz = s->pz[i] - cz, prad = s->pr[i];
        double relc = tx*v[0] + ty*v[1] + tz*v[2], cr;
        if (fabs(relc) < 1e-9) cr = prad;
        else if (fabs(relc) < prad && fabs(relc) < 0.9*sample) cr = sqrt(prad*prad - relc*relc);
        else { cr = -1e10; act[i] = 0; }
        circrad[i] = cr;
        if (cr > endrad) { s->requiv = 1.0e6; s->escaped = 1; goto done; }
    }
    double emax = -1e10, emin = 1e10, nmax = -1e10, nmin = 1e10;
    int *vi = malloc((n > 0 ? n : 1) * sizeof(int)); int nv = 0;
    for (int i = 0; i < n; i++) {
        if (!act[i] || circrad[i] >= endrad) continue;
        double tx = s->px[i] - cx, ty = s->py[i] - cy, tz = s->pz[i] - cz;
        double ee = tx*e[0] + ty*e[1] + tz*e[2], q = tx*nn[0] + ty*nn[1] + tz*nn[2];
        ec[i] = ee; nc[i] = q;
        double cr = circrad[i];
        if (ee + cr > emax) emax = ee + cr;
        if (ee - cr < emin) emin = ee - cr;
        if (q + cr > nmax) nmax = q + cr;
        if (q - cr < nmin) nmin = q - cr;
        vi[nv++] = i;
    }
    if (fabs(emax + 1e10) < 0.001) { s->requiv = 0.0; free(vi); goto done; }
    double size1 = 0.25, area = 0.0;
    int tripE = (int)((emax - emin + size1) / size1); if (tripE < 0) tripE = 0;
    int tripN = (int)((nmax - nmin + size1) / size1); if (tripN < 0) tripN = 0;
    /* cell origins exactly as the original forms them: repeated addition */
    double *edc = malloc((tripE + 1) * sizeof(double)), *ndc = malloc((tripN + 1) * sizeof(double));
    { double a = emin; for (int ie = 0; ie <= tripE; ie++) { edc[ie] = a; a += size1; } }
    { double a = nmin; for (int jn = 0; jn <= tripN; jn++) { ndc[jn] = a; a += size1; } }
    long ncells = (long)tripE * tripN;
    ho_area A = { nv, vi, ec, nc, circrad, edc, ndc, size1, 0.5*root2*size1,
                  calloc(ncells > 0 ? ncells : 1, 1), tripN, NULL,
                  malloc((ncells > 0 ? ncells : 1) * sizeof(int)),
                  malloc((ncells > 0 ? ncells : 1) * sizeof(int)), 0, 0 };
    if (tripE > 0 && tripN > 0) ho_block(&A, 0, tripE, 0, tripN, vi, nv);
    ho_cell *curl = malloc((HO_SNMAX + 8) * sizeof(ho_cell)), *newl = malloc((HO_SNMAX + 8) * sizeof(ho_cell));
    int ncur = 0;
    for (int ie = 0; ie < tripE; ie++)
        for (int jn = 0; jn < tripN; jn++) {
            int cell = ie * tripN + jn;
            if (A.mask[cell] == 2) area += size1*size1;
            else if (A.mask[cell] == 1 && ncur < HO_SNMAX + 8) {
                curl[ncur].bx = edc[ie]; curl[ncur].by = ndc[jn]; curl[ncur].root = cell; ncur++;
            }
        }
    double areag = (double)ncur * size1 * size1;
    double requiv = sqrt((area + 0.5*areag) / pi);
    double cur = size1;
    int ncycle = 0;
    for (;;) {
        ncycle++;
        double nsz = 0.5 * cur; int newCount = 0, nnew = 0;
        for (int c = 0; c < ncur; c++) {
            double bx = curl[c].bx, by = curl[c].by;
            const int *cl = A.cand + A.cand_off[curl[c].root]; int ncl = A.cand_n[curl[c].root];
            double edgeE[8] = {bx, bx+0.5*cur, bx+cur, bx, bx+cur, bx, bx+0.5*cur, bx+cur};
            double edgeN[8] = {by, by, by, by+0.5*cur, by+0.5*cur, by+cur, by+cur, by+cur};
            int edblk[8] = {0,0,0,0,0,0,0,0};
            double cenE = bx + 0.5*cur, cenN = by + 0.5*cur;
            int sumall = 0;
            for (int q = 0; q < ncl; q++) {
                int i = cl[q];
                double cr = circrad[i], de = cenE - ec[i], dn = cenN - nc[i];
                double dist2 = de*de + dn*dn, w = cr + 0.5*root2*cur;
                if (dist2 > w*w) continue;
                for (int ei = 0; ei < 8; ei++) {
                    double de2 = edgeE[ei] - ec[i], dn2 = edgeN[ei] - nc[i];
                    if (de2*de2 + dn2*dn2 < cr*cr) edblk[ei] = 1;
                }
                sumall = 0; for (int ei = 0; ei < 8; ei++) sumall += edblk[ei];
                if (sumall == 8) break;
            }
            sumall = 0; for (int ei = 0; ei < 8; ei++) sumall += edblk[ei];
            if (sumall == 8) area += cur*cur;
            else if (sumall == 0) { }
            else {
                if (newCount + 4 >= HO_SNMAX) newCount += 4;
                else {
                    newCount += 4;
                    int root = curl[c].root;
                    newl[nnew].bx = bx;       newl[nnew].by = by;       newl[nnew].root = root; nnew++;
                    newl[nnew].bx = bx + nsz; newl[nnew].by = by;       newl[nnew].root = root; nnew++;
                    newl[nnew].bx = bx;       newl[nnew].by = by + nsz; newl[nnew].root = root; nnew++;
                    newl[nnew].bx = bx + nsz; newl[nnew].by = by + nsz; newl[nnew].root = root; nnew++;
                }
            }
        }
        areag = (double)newCount * nsz * nsz;
        double oldreq = requiv;
        requiv = sqrt((area + 0.5*areag) / pi);
        if (fabs(requiv - oldreq) > 0.0005 || ncycle < 4) {
            if (newCount + 4 <= HO_SNMAX) {
                cur = nsz;
                ho_cell *t = curl; curl = newl; newl = t; ncur = nnew;
                continue;
            }
            break;
        }
        break;
    }
    s->requiv = requiv;
    free(curl); free(newl); free(A.mask); free(A.cand); free(A.cand_off); free(A.cand_n);
    free(edc); free(ndc); free(vi);
done:
    free(circrad); free(ec); free(nc);
}

/* addend.f: escape spheres beyond the second-to-last slice at each end. */
static int ho_addend(double lcx, double lcy, double lcz, double sample_signed,
                     const double v[3], const double e[3], const double nn[3],
                     double endrad, double *ox, double *oy, double *oz, double *orad) {
    double ssign = sample_signed > 0.0 ? 1.0 : -1.0, half = 0.5 * endrad;
    int k = 0;
    for (int zc = 0; zc <= 4; zc++) {
        double zd = ssign * (double)zc * half;
        for (int ncount = -2; ncount <= 2; ncount++) {
            double nd = (double)ncount * half;
            for (int ecount = -2; ecount <= 2; ecount++) {
                double ed = (double)ecount * half;
                double tx = lcx + zd*v[0] + nd*nn[0] + ed*e[0];
                double ty = lcy + zd*v[1] + nd*nn[1] + ed*e[1];
                double tz = lcz + zd*v[2] + nd*nn[2] + ed*e[2];
                double c = clearance(tx, ty, tz);
                if (c > endrad) { ox[k] = tx; oy[k] = ty; oz[k] = tz; orad[k] = c; k++; }
            }
        }
    }
    return k;
}

/* ------------------------------------------------------------- emitters */

/* Slice bookkeeping: the march stores slices in ascending t. HOLE numbers
   the start slice 0, +ve slices 1.., -ve slices -1.. downward. */
typedef struct {
    int n;            /* slices */
    int *irec;        /* HOLE record number per slice */
    int *order;       /* slice indices in HOLE write order (pos asc, then neg desc) */
    int npos, nneg;   /* counts */
    double sample;
    double v[3], e[3], nn[3];
    ho_slice *conn;   /* per-slice Connolly result or NULL */
    double probe;
} ho_layout;

static void ho_layout_init(ho_layout *L, double sample) {
    L->n = nslice; L->sample = sample; L->conn = NULL; L->probe = 0;
    for (int i = 0; i < 3; i++) L->v[i] = W[i];
    ho_calper(L->v, L->e, L->nn);
    L->irec = malloc(nslice * sizeof(int));
    L->order = malloc(nslice * sizeof(int));
    L->npos = L->nneg = 0;
    for (int i = 0; i < nslice; i++) {
        L->irec[i] = (int)lround(slt[i] / sample);
        if (slt[i] >= -1e-6) L->npos++; else L->nneg++;
    }
    int k = 0;
    for (int i = 0; i < nslice; i++) if (slt[i] >= -1e-6) L->order[k++] = i;
    for (int i = nslice - 1; i >= 0; i--) if (slt[i] < -1e-6) L->order[k++] = i;
}

static void ho_layout_free(ho_layout *L) {
    free(L->irec); free(L->order);
    if (L->conn) { for (int i = 0; i < L->n; i++) ho_slice_free(&L->conn[i]); free(L->conn); }
}

/* Connolly pass over every slice, parallel. */
static void ho_run_conn(ho_layout *L, double endrad, double probe, double grid) {
    L->conn = calloc(L->n, sizeof(ho_slice));
    L->probe = probe;
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < L->n; i++) {
        ho_slice *s = &L->conn[i];
        if (!ho_concal(slx[i], sly[i], slz[i], L->e, L->nn, endrad, probe, grid, s)) continue;
        ho_coarea(s, slx[i], sly[i], slz[i], L->e, L->nn, L->v, L->sample, endrad);
    }
}

static void ho_write_slice(FILE *f, const ho_layout *L, int i, double endrad) {
    if (!L->conn) { ho_sph_line(f, L->irec[i], slx[i], sly[i], slz[i], slr[i], slr[i]); return; }
    const ho_slice *s = &L->conn[i];
    if (s->fallback) return;
    double wreq = s->requiv > 1000.0 ? 999.99 : s->requiv;
    int wresno = L->irec[i];
    for (int k = 0; k < s->n; k++) {
        if (s->act[k] || s->pr[k] > endrad) {
            ho_sph_line(f, wresno, s->px[k], s->py[k], s->pz[k], s->pr[k], wreq);
            if (s->pr[k] > endrad) fputs("LAST-REC-END\n", f);
            wresno = -999;
        }
    }
}

static void ho_write_addend(FILE *f, const ho_layout *L, int from, double ssample, double endrad) {
    double ox[125], oy[125], oz[125], orad[125];
    int k = ho_addend(slx[from], sly[from], slz[from], ssample, L->v, L->e, L->nn, endrad, ox, oy, oz, orad);
    for (int j = 0; j < k; j++) { ho_sph_line(f, -888, ox[j], oy[j], oz[j], orad[j], 0.0); fputs("LAST-REC-END\n", f); }
}

static int ho_write_sph(const char *path, const ho_layout *L, double endrad) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return 0; }
    if (L->npos == 0) { fclose(f); return 0; }
    for (int k = 0; k < L->npos; k++) ho_write_slice(f, L, L->order[k], endrad);
    if (L->npos >= 2) ho_write_addend(f, L, L->order[L->npos - 2], L->sample, endrad);
    int i0 = L->order[0];
    ho_sph_line(f, L->irec[i0], slx[i0], sly[i0], slz[i0], slr[i0], slr[i0]);
    for (int k = L->npos; k < L->n; k++) ho_write_slice(f, L, L->order[k], endrad);
    if (L->nneg >= 2) ho_write_addend(f, L, L->order[L->n - 2], -L->sample, endrad);
    fclose(f);
    return 1;
}

/* Profile TSV: slice rows and mid-point rows sorted by coord. */
/* HOLE's profile stops at ENDRAD: the slice that crossed it is kept in the
   .sph (it becomes a mouth clip sphere) but never printed as a profile row.
   Printing it made the r^2-weighted Volume readout differ by ~20% between the
   two search engines on a pore whose bottleneck agreed to 0.003%. */
static int ho_write_tsv(const char *path, const ho_layout *L, double endrad) {
    int n = L->n, m = 2 * n - 1;
    if (n < 1) return 0;
    double *coord = malloc(m * sizeof(double)), *rad = malloc(m * sizeof(double)),
           *x = malloc(m * sizeof(double)), *y = malloc(m * sizeof(double)), *z = malloc(m * sizeof(double));
    int *ri = malloc(m * sizeof(int));
    int k = 0;
    for (int i = 0; i < n; i++) {
        if (i > 0) {   /* slices are ascending in t already, so neighbours are adjacent */
            double mx = 0.5*(slx[i-1] + slx[i]), my = 0.5*(sly[i-1] + sly[i]), mz = 0.5*(slz[i-1] + slz[i]);
            x[k] = mx; y[k] = my; z[k] = mz;
            coord[k] = mx*L->v[0] + my*L->v[1] + mz*L->v[2];
            rad[k] = clearance(mx, my, mz); ri[k] = -1; k++;
        }
        x[k] = slx[i]; y[k] = sly[i]; z[k] = slz[i];
        coord[k] = slx[i]*L->v[0] + sly[i]*L->v[1] + slz[i]*L->v[2];
        rad[k] = slr[i]; ri[k] = i; k++;
    }
    /* coord of the start slice, arc length zero there */
    double t0 = 0; for (int i = 0; i < n; i++) if (L->irec[i] == 0) t0 = slx[i]*L->v[0] + sly[i]*L->v[1] + slz[i]*L->v[2];
    double cum = 0, zero = 0, zbest = 1e30;
    double *dist = malloc(m * sizeof(double));
    for (int j = 0; j < m; j++) {
        if (j > 0) cum += sqrt((x[j]-x[j-1])*(x[j]-x[j-1]) + (y[j]-y[j-1])*(y[j]-y[j-1]) + (z[j]-z[j-1])*(z[j]-z[j-1]));
        dist[j] = cum;
        if (fabs(coord[j] - t0) < zbest) { zbest = fabs(coord[j] - t0); zero = cum; }
    }
    const double pi = 2.0 * asin(1.0), samp = fabs(L->sample);
    double cratio = 0;
    if (L->conn) {
        double cs = 0; int cn = 0;
        for (int i = 0; i < n; i++) {
            double rq = L->conn[i].requiv;
            if (rq < 1.0e6 && slr[i] != 0) { cs += rq / slr[i]; cn++; }
        }
        if (cn > 0) cratio = cs / cn;
    }
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return 0; }
    fputs("coord\tradius\tcen_line_d\tsum_s_over_area\trequiv\tconn_s_over_area\trequiv_estim\tcap_rad\n", f);
    double F = 0, CF = 0;
    for (int j = 0; j < m; j++) {
        if (rad[j] > endrad) continue;
        if (fabs(rad[j]) > 1e-9) F += 0.5 * samp / (pi * rad[j] * rad[j]);
        char c5[32] = "", c6[32] = "", c7[32] = "";
        if (L->conn && ri[j] >= 0) {
            double rq = L->conn[ri[j]].requiv;
            double re = rq > 0.999e6 ? cratio * rad[j] : rq;
            if (re > 0.0) CF += samp / (pi * re * re);
            snprintf(c5, sizeof c5, "%.3f", rq);
            snprintf(c6, sizeof c6, "%.3f", CF);
            snprintf(c7, sizeof c7, "%.3f", re);
        }
        fprintf(f, "%.5f\t%.5f\t%.5f\t%.5f\t%s\t%s\t%s\t\n", coord[j], rad[j], dist[j] - zero, F, c5, c6, c7);
    }
    fclose(f);
    free(coord); free(rad); free(x); free(y); free(z); free(ri); free(dist);
    return 1;
}

/* ------------------------------------------------------- ignore + centres */

/* HOLE's IGNORE card: residue names (3 chars, case-insensitive) dropped
   before the search. Applied to the loaded atom arrays in place. */
static int ho_apply_ignore(const char *spec, const char (*resnames)[4]) {
    char names[64][4]; int nn = 0;
    const char *p = spec;
    while (*p && nn < 64) {
        while (*p == ',' || *p == ' ') p++;
        if (!*p) break;
        int k = 0;
        while (*p && *p != ',' && *p != ' ') { if (k < 3) names[nn][k++] = (char)toupper((unsigned char)*p); p++; }
        names[nn][k] = 0; nn++;
    }
    if (nn == 0) return 0;
    int kept = 0, dropped = 0;
    for (int i = 0; i < natom; i++) {
        char r[4] = {0,0,0,0}; int k = 0;
        for (int j = 0; j < 3; j++) if (resnames[i][j] != ' ' && resnames[i][j]) r[k++] = (char)toupper((unsigned char)resnames[i][j]);
        int drop = 0;
        for (int q = 0; q < nn; q++) if (!strcmp(r, names[q])) { drop = 1; break; }
        if (drop) { dropped++; continue; }
        ax[kept] = ax[i]; ay[kept] = ay[i]; az[kept] = az[i]; ar[kept] = ar[i]; kept++;
    }
    natom = kept;
    return dropped;
}

/* Validation aid: take the centreline from a file of "t x y z r" lines
   (t relative to CPOINT) instead of searching, so the Connolly pass can be
   checked against HOLE's own dots on HOLE's own centres. */
static int ho_load_centres(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return 0; }
    char line[256]; nslice = 0;
    while (fgets(line, sizeof line, f) && nslice < MAXSL) {
        double t, x, y, z, r;
        if (sscanf(line, "%lf %lf %lf %lf %lf", &t, &x, &y, &z, &r) != 5) continue;
        slt[nslice] = t; slx[nslice] = x; sly[nslice] = y; slz[nslice] = z; slr[nslice] = r; nslice++;
    }
    fclose(f);
    /* ascending t */
    for (int i = 1; i < nslice; i++) {
        double kt = slt[i], kx = slx[i], ky = sly[i], kz = slz[i], kr = slr[i]; int j = i - 1;
        while (j >= 0 && slt[j] > kt) { slt[j+1] = slt[j]; slx[j+1] = slx[j]; sly[j+1] = sly[j]; slz[j+1] = slz[j]; slr[j+1] = slr[j]; j--; }
        slt[j+1] = kt; slx[j+1] = kx; sly[j+1] = ky; slz[j+1] = kz; slr[j+1] = kr;
    }
    return nslice;
}
