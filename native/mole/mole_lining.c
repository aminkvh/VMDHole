/* MOLE 2's lining layers and physico-chemical properties.
 *
 * Tunnel.GetLiningLayers samples the centreline three times per Angstrom, takes
 * the five nearest atoms at each sample, groups them into residues, then merges
 * consecutive samples whose residue set is identical into one layer carrying the
 * MINIMUM of the widths over its run. Tunnel.ComputePhysicoChemicalProperties
 * then reduces those layers to per-tunnel numbers.
 *
 * Faithfulness notes, because several of these read as bugs:
 *  - a residue counts as backbone-touched only if EVERY one of its atoms among
 *    the five is a backbone atom (Tunnel.cs:228);
 *  - the backbone half of the property sum uses ASN's polarity and GLY's
 *    hydrophobicity/hydropathy, not BACKBONE's own (PhysicoChemicalProperties.cs
 *    :238-247) - BACKBONE supplies only LogP/LogD/LogS;
 *  - that same loop increments the denominator unconditionally while the
 *    side-chain loop increments it only on a table hit;
 *  - layers are ordered by residue NUMBER alone, with no chain tiebreak.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_lining.h"
#include "../xalloc.h"

/* TunnelPhysicoChemicalPropertyTable. Charge, Ionizable, Hydropathy,
   Hydrophobicity, Polarity, LogP, LogD, LogS, Mutability. */
typedef struct {
    const char *resn;
    int    charge, ionizable, mutability;
    double hydropathy, hydrophobicity, polarity, logp, logd, logs;
} res_info;

static const res_info TABLE[] = {
    { "ALA",  0, 0, 100,  1.8,  0.02,  0.00,  1.08,  1.08,  0.59 },
    { "ARG",  1, 1,  83, -4.5, -0.42, 52.00, -0.08, -2.49,  1.63 },
    { "ASN",  0, 0, 104, -3.5, -0.77,  3.38, -1.03, -1.03,  0.54 },
    { "ASP", -1, 1,  86, -3.5, -1.04, 49.70, -0.22, -3.00,  2.63 },
    { "CYS",  0, 0,  44,  2.5,  0.77,  1.48,  0.84,  0.84,  0.16 },
    { "GLU", -1, 1,  77, -3.5, -1.14, 49.90,  0.48, -2.12,  2.23 },
    { "GLN",  0, 0,  84, -3.5, -1.10,  3.53, -0.33, -0.33,  0.13 },
    { "GLY",  0, 0,  50, -0.4, -0.80,  0.00,  0.00,  0.00,  0.00 },
    { "HIS",  0, 0,  91, -3.2,  0.26, 51.60, -0.01, -0.11, -0.20 },
    { "ILE",  0, 0, 103,  4.5,  1.81,  0.13,  2.24,  2.24, -1.85 },
    { "LEU",  0, 0,  54,  3.8,  1.14,  0.13,  2.08,  2.08, -1.79 },
    { "LYS",  1, 1,  72, -3.9, -0.41, 49.50,  0.70, -1.91,  1.46 },
    { "MET",  0, 0,  93,  1.9,  1.00,  1.43,  1.48,  1.48, -0.72 },
    { "PHE",  0, 0,  51,  2.8,  1.35,  0.35,  2.49,  2.49, -1.81 },
    { "PRO",  0, 0,  58, -1.6, -0.09,  1.58,  1.80,  1.80, -1.30 },
    { "SER",  0, 0, 117, -0.8, -0.97,  1.67, -0.52, -0.52,  1.11 },
    { "THR",  0, 0, 107, -0.7, -0.77,  1.66, -0.16, -0.16,  0.77 },
    { "TRP",  0, 0,  25, -0.9,  1.71,  2.10,  2.59,  2.59, -2.48 },
    { "TYR",  0, 0,  50, -1.3,  1.11,  1.61,  2.18,  2.18, -1.44 },
    { "VAL",  0, 0,  98,  4.2,  1.13,  0.13,  1.80,  1.80, -1.30 }
};

/* The constants the backbone half of the sum contributes. Deliberately three
   different table rows - see the header comment. */
#define BB_POLARITY        3.38     /* ASN */
#define BB_HYDROPHOBICITY (-0.80)   /* GLY */
#define BB_HYDROPATHY     (-0.40)   /* GLY */
#define BB_LOGP           (-0.86)   /* BACKBONE */
#define BB_LOGD           (-0.86)   /* BACKBONE */
#define BB_LOGS             0.81    /* BACKBONE */

static const res_info *lookup(const char *resn)
{
    size_t i, j;
    char u[8];
    for (i = 0; i + 1 < sizeof u && resn[i]; i++)
        u[i] = (resn[i] >= 'a' && resn[i] <= 'z') ? (char)(resn[i] - 32) : resn[i];
    u[i] = 0;
    for (j = 0; j < sizeof TABLE / sizeof TABLE[0]; j++)
        if (!strcmp(u, TABLE[j].resn)) return &TABLE[j];
    return NULL;    /* GetResidueProperties returns null and the caller skips */
}

int mole_pivot_residues(const mole_atoms *A, int *pres, mole_residues *R)
{
    int i, k, p, np = 0;
    R->n = 0;
    /* PIVOT order, not file order - these are per-pivot arrays. Requires
       mole_pivots to have run. */
    for (p = 0; p < A->npiv; p++) {
        i = A->order[p];
        /* A residue's atoms are contiguous in pivot order too, so the last one
           entered is almost always the answer. Without this first probe the
           scan below is quadratic in the atom count. */
        k = R->n - 1;
        if (k < 0 || R->seq[k] != A->seq[i] || strcmp(R->chain[k], A->chain[i])
            || strcmp(R->name[k], A->resn[i])) {
            for (k = 0; k < R->n; k++)
                if (R->seq[k] == A->seq[i] && !strcmp(R->chain[k], A->chain[i])
                    && !strcmp(R->name[k], A->resn[i])) break;
        }
        if (k == R->n) {
            /* REFUSE rather than truncate, matching mole_read_atoms' policy at
               the atom cap: silently dropping residues does not fail, it
               produces a plausible-but-wrong answer - layers lose residues from
               their averaged charge/hydropathy/logP rows, cavity boundary sets
               shrink, and an emptied boundary set relabels a Cavity as a Void. */
            if (R->n >= MOLE_MAXRES) {
                fprintf(stderr, "mole_pivot_residues: more than %d residues - "
                        "refusing to analyse a truncated structure. Narrow the "
                        "selection or rebuild with a larger MOLE_MAXRES.\n",
                        MOLE_MAXRES);
                return -1;
            }
            strcpy(R->name[k], A->resn[i]);
            strcpy(R->chain[k], A->chain[i]);
            R->seq[k] = A->seq[i];
            R->n++;
        }
        pres[np++] = k;
    }
    return np;
}

/* ---- one sample's residue set ------------------------------------------ */

typedef struct {
    int  nres;
    int  res[MOLE_LAYER_MAXRES];
    char isbb[MOLE_LAYER_MAXRES];
} sample_lining;

/* Group the five nearest atoms by residue, then order by residue NUMBER.
   MOLE's OrderBy is stable and has no chain tiebreak, so equal numbers keep
   the order the k-d tree returned them in. */
static void group_residues(const int *sel, int nsel, const int *pres,
                           const mole_residues *R, const int *pbb,
                           sample_lining *out)
{
    int i, k;
    out->nres = 0;
    for (i = 0; i < nsel; i++) {
        int r = pres[sel[i]];
        if (r < 0) continue;
        for (k = 0; k < out->nres; k++) if (out->res[k] == r) break;
        if (k == out->nres) { out->res[k] = r; out->isbb[k] = 1; out->nres++; }
        /* IsBackbone = g.All(a => a.IsBackboneAtom()) - one side-chain atom
           among the five is enough to make the whole residue non-backbone. */
        if (!pbb[sel[i]]) out->isbb[k] = 0;
    }
    for (i = 1; i < out->nres; i++) {          /* stable insertion sort */
        int rv = out->res[i]; char bv = out->isbb[i];
        for (k = i; k > 0 && R->seq[out->res[k-1]] > R->seq[rv]; k--) {
            out->res[k] = out->res[k-1]; out->isbb[k] = out->isbb[k-1];
        }
        out->res[k] = rv; out->isbb[k] = bv;
    }
}

/* Collect one half of a sample's residues, in order. */
static int half(const sample_lining *s, int wantbb, int *out)
{
    int i, n = 0;
    for (i = 0; i < s->nres; i++) if (s->isbb[i] == wantbb) out[n++] = s->res[i];
    return n;
}

static int same_lining(const sample_lining *a, const sample_lining *b)
{
    /* MOLE keys a layer on the non-backbone identifiers concatenated, then "B",
       then the backbone ones - so the SPLIT matters as much as the set, and the
       comparison is ordered. Element-wise here rather than by concatenated
       string; the two differ only if two identifier lists concatenate alike. */
    int i, j, ka[MOLE_LAYER_MAXRES], kb[MOLE_LAYER_MAXRES], na, nb;
    for (j = 0; j < 2; j++) {
        na = half(a, j, ka);
        nb = half(b, j, kb);
        if (na != nb) return 0;
        for (i = 0; i < na; i++) if (ka[i] != kb[i]) return 0;
    }
    return 1;
}

/* ---- properties -------------------------------------------------------- */

/* CalculateHydrophibilicyPolarityHydropathy over a residue list already split
   into side-chain and backbone-touched. nbb residues contribute the constants
   AND the denominator with no table lookup at all. */
static void hydro_block(const int *res, const char *isbb, int nres,
                        const mole_residues *R, mole_props *P)
{
    int i, count = 0;
    P->hydropathy = P->hydrophobicity = P->polarity = 0.0;
    P->logp = P->logd = P->logs = 0.0;
    for (i = 0; i < nres; i++) {
        const res_info *info;
        if (isbb[i]) continue;
        info = lookup(R->name[res[i]]);
        if (!info) continue;
        count++;
        P->hydropathy += info->hydropathy;
        P->hydrophobicity += info->hydrophobicity;
        P->polarity += info->polarity;
        P->logp += info->logp; P->logd += info->logd; P->logs += info->logs;
    }
    for (i = 0; i < nres; i++) {
        if (!isbb[i]) continue;
        count++;
        P->polarity += BB_POLARITY;
        P->hydrophobicity += BB_HYDROPHOBICITY;
        P->hydropathy += BB_HYDROPATHY;
        P->logp += BB_LOGP; P->logd += BB_LOGD; P->logs += BB_LOGS;
    }
    if (count == 0) {
        P->hydropathy = P->hydrophobicity = P->polarity = 0.0;
        P->logp = P->logd = P->logs = 0.0;
    } else {
        P->hydropathy /= count; P->hydrophobicity /= count;
        P->polarity /= count;
        P->logp /= count; P->logd /= count; P->logs /= count;
    }
}

/* Charge, Ionizable and Mutability come from the side-chain residues only. */
static void count_block(const int *res, const char *isbb, int nres,
                        const mole_residues *R, mole_props *P)
{
    int i, count = 0;
    double mut = 0.0;
    P->charge = P->ionizable = P->npos = P->nneg = 0;
    for (i = 0; i < nres; i++) {
        const res_info *info;
        if (isbb[i]) continue;
        info = lookup(R->name[res[i]]);
        if (!info) continue;
        count++;
        P->charge += info->charge;
        P->ionizable += info->ionizable;
        if (info->charge > 0) P->npos++;
        else if (info->charge < 0) P->nneg++;
        mut += info->mutability;
    }
    if (count) mut /= count;
    P->mutability = (int)mut;
}

/* ---- ResidueFlow ------------------------------------------------------- */

/* TunnelLining.ComputeFlow. Walks the layers in order and, WITHIN each layer,
   the residues sorted by chain then number - which is not the order the layer
   itself keeps (that one has no chain tiebreak). The key carries the backbone
   flag, so one residue can occupy two flow slots. */
int mole_build_flow(mole_layer *L, int nl, const mole_residues *R,
                    mole_flow_entry *flow)
{
    int i, j, k, n = 0;
    for (i = 0; i < nl; i++) {
        int ord[MOLE_LAYER_MAXRES], nr = L[i].nres;
        for (j = 0; j < nr; j++) ord[j] = j;
        for (j = 1; j < nr; j++) {          /* stable, by (chain, number) */
            int v = ord[j];
            for (k = j; k > 0; k--) {
                int a = ord[k-1], c = mole_chain_cmp(R->chain[L[i].res[a]],
                                                    R->chain[L[i].res[v]]);
                if (c < 0 || (c == 0 && R->seq[L[i].res[a]] <= R->seq[L[i].res[v]])) break;
                ord[k] = ord[k-1];
            }
            ord[k] = v;
        }
        for (j = 0; j < nr; j++) {
            int s = ord[j];
            for (k = 0; k < n; k++)
                if (flow[k].res == L[i].res[s] && flow[k].isbb == L[i].isbb[s]) break;
            if (k == n) { flow[n].res = L[i].res[s]; flow[n].isbb = L[i].isbb[s]; n++; }
            L[i].flow[s] = k;
        }
    }
    return n;
}

/* ---- FindBottleneck / UpdateMinima ------------------------------------- */

static void mark_minima(mole_layer *L, int nl)
{
    int i, first = -1, best = -1;
    if (nl <= 0) return;
    /* SkipWhile(EndDistance < 3), then the FIRST layer achieving the minimum
       radius. With no such layer MOLE takes the min over layers[1..] and,
       unlike the normal path, does NOT mark it. */
    for (i = 0; i < nl; i++) if (L[i].end >= 3.0) { first = i; break; }
    if (first < 0 && getenv("MOLE_BRANCH_AUDIT"))
        fprintf(stderr, "BRANCH findbottleneck-no-candidates nl=%d\n", nl);
    if (first >= 0) {
        best = first;
        for (i = first; i < nl; i++) if (L[i].radius < L[best].radius) best = i;
        L[best].localmin = 1;
    }
    for (i = 1; i < nl - 1; i++)
        if (L[i].radius < L[i-1].radius && L[i].radius < L[i+1].radius)
            L[i].localmin = 1;
    /* Runs on the values just written, so it propagates from the bottleneck too. */
    for (i = 1; i < nl; i++)
        if (L[i-1].localmin && fabs(L[i].radius - L[i-1].radius) < 0.005)
            L[i].localmin = 1;
    if (nl > 1 && L[nl-1].radius < L[nl-2].radius) L[nl-1].localmin = 1;
}

/* ---- build ------------------------------------------------------------- */

int mole_lining_build(const mole_tunnel_profile *p,
                      const double *axyz, const double *arad, int na,
                      const int *pres, const mole_residues *R, const int *pbb,
                      mole_lining *out)
{
    /* samplesPerAngstrom = 3, numSurroundingAtoms = 5 (Tunnel.cs:194). */
    int n = (int)ceil(p->length * 3.0), i, k, nl = 0, start = 0;
    double dt, px = 0, py = 0, pz = 0, dist = 0;
    sample_lining *S = NULL, cur;
    double *sd = NULL, *sr = NULL, *sf = NULL, *sb = NULL;
    mole_layer *L = NULL;
    double wsum = 0.0;
    mole_props w;
    int *allres = NULL; char *allbb = NULL; int nall = 0;

    memset(out, 0, sizeof *out);
    if (n < 2) n = 2;
    dt = 1.0 / (n - 1);
    S  = xa_malloc((size_t)n * sizeof *S);
    sd = xa_malloc((size_t)n * sizeof *sd);
    sr = xa_malloc((size_t)n * sizeof *sr);
    sf = xa_malloc((size_t)n * sizeof *sf);
    sb = xa_malloc((size_t)n * sizeof *sb);
    L  = xa_malloc((size_t)n * sizeof *L);
    /* calloc, not malloc: nall is bounded by n*MOLE_LAYER_MAXRES but nothing in
       the code says so, and a half-filled buffer read past the fill is exactly
       the class of bug that is invisible until a structure triggers it. */
    allres = calloc((size_t)n * MOLE_LAYER_MAXRES, sizeof *allres);
    allbb  = calloc((size_t)n * MOLE_LAYER_MAXRES, 1);
    if (!S || !sd || !sr || !sf || !sb || !L || !allres || !allbb) goto fail;

    for (i = 0; i < n; i++) {
        double t = dt * i, x, y, z;
        int sel[5], nsel;
        x = mole_spline_eval(&p->sx, t);
        y = mole_spline_eval(&p->sy, t);
        z = mole_spline_eval(&p->sz, t);
        if (i > 0) {
            double ddx = x - px, ddy = y - py, ddz = z - pz;
            dist += sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
        }
        /* The t=0 centreline point at full precision. It is the first control
           path tetrahedron's circumcentre exactly - the spline returns its
           first knot - and it is where equidistant atoms are ordered by
           round-off. */
        if (i == 0 && getenv("MOLE_T0_DUMP")) {
            FILE *tf = fopen(getenv("MOLE_T0_DUMP"), "a");
            if (tf) { fprintf(tf, "%d %.17g %.17g %.17g\n", n, x, y, z); fclose(tf); }
        }
        px = x; py = y; pz = z;
        { double q[3]; q[0] = x; q[1] = y; q[2] = z;
          nsel = mole_nearest5(q, axyz, arad, na, sel); }
        group_residues(sel, nsel, pres, R, pbb, &S[i]);
        sd[i] = dist;
        sr[i] = mole_spline_eval(&p->sr, t);
        sf[i] = mole_spline_eval(&p->sfr, t);
        sb[i] = mole_spline_eval(&p->sbr, t);
    }

    if (getenv("MOLE_LAYER_DUMP")) {
        FILE *lf = fopen(getenv("MOLE_LAYER_DUMP"), "a");
        if (lf) {
            fprintf(lf, "TUNNEL samples=%d\n", n);
            for (i = 0; i < n; i++) {
                int q, firstb = 1;
                fprintf(lf, "  %d nb=[", i);
                for (q = 0; q < S[i].nres; q++) {
                    if (S[i].isbb[q]) continue;
                    fprintf(lf, "%s%d %s", firstb ? "" : ",",
                            R->seq[S[i].res[q]], R->chain[S[i].res[q]]);
                    firstb = 0;
                }
                fprintf(lf, "] bb=[");
                firstb = 1;
                for (q = 0; q < S[i].nres; q++) {
                    if (!S[i].isbb[q]) continue;
                    fprintf(lf, "%s%d %s", firstb ? "" : ",",
                            R->seq[S[i].res[q]], R->chain[S[i].res[q]]);
                    firstb = 0;
                }
                fprintf(lf, "]\n");
            }
            fclose(lf);
        }
    }

    /* Merge equal-lining runs. Each emitted layer keeps the FIRST sample's
       residues and distance but the MINIMUM width over its whole run. */
    cur = S[0];
    for (i = 1; i <= n; i++) {
        if (i < n && same_lining(&S[i], &cur)) continue;
        L[nl].nres = cur.nres;
        for (k = 0; k < cur.nres; k++) { L[nl].res[k] = cur.res[k]; L[nl].isbb[k] = cur.isbb[k]; }
        L[nl].start = sd[start];
        L[nl].end   = (i < n) ? sd[i] : sd[n-1];
        L[nl].radius = sr[start]; L[nl].freeradius = sf[start]; L[nl].bradius = sb[start];
        for (k = start; k < i; k++) {
            if (sr[k] < L[nl].radius) L[nl].radius = sr[k];
            if (sf[k] < L[nl].freeradius) L[nl].freeradius = sf[k];
            if (sb[k] < L[nl].bradius) L[nl].bradius = sb[k];
        }
        L[nl].localmin = 0;
        count_block(L[nl].res, L[nl].isbb, L[nl].nres, R, &L[nl].props);
        hydro_block(L[nl].res, L[nl].isbb, L[nl].nres, R, &L[nl].props);
        nl++;
        if (i < n) { cur = S[i]; start = i; }
    }
    out->flow = xa_malloc((size_t)nl * MOLE_LAYER_MAXRES * sizeof *out->flow);
    if (!out->flow) goto fail;
    out->nflow = mole_build_flow(L, nl, R, out->flow);
    mark_minima(L, nl);

    /* Whole tunnel. The counting half deduplicates residues across layers; the
       hydro half deliberately does not, so a residue lining ten layers is
       weighted ten times there and once here. */
    for (i = 0; i < nl; i++)
        for (k = 0; k < L[i].nres; k++) {
            allres[nall] = L[i].res[k]; allbb[nall] = L[i].isbb[k]; nall++;
        }
    hydro_block(allres, allbb, nall, R, &out->props);
    {
        int m = 0, j;
        for (i = 0; i < nall; i++) {
            if (allbb[i]) continue;
            for (j = 0; j < m; j++) if (allres[j] == allres[i]) break;
            if (j == m) { allres[m] = allres[i]; allbb[m] = 0; m++; }
        }
        count_block(allres, allbb, m, R, &out->props);
    }

    /* LayerWeightedPhysicoChemicalProperties: w = EndDistance - StartDistance.
       Charge/Ionizable/counts are carried over unweighted. */
    memset(&w, 0, sizeof w);
    {
        double wm = 0.0;    /* Mutability is an int per layer, summed as a double */
        for (i = 0; i < nl; i++) {
            double ww = L[i].end - L[i].start;
            wsum += ww;
            w.hydropathy += ww * L[i].props.hydropathy;
            w.hydrophobicity += ww * L[i].props.hydrophobicity;
            w.polarity += ww * L[i].props.polarity;
            w.logp += ww * L[i].props.logp;
            w.logd += ww * L[i].props.logd;
            w.logs += ww * L[i].props.logs;
            wm += ww * L[i].props.mutability;
        }
        /* Below the guard MOLE leaves the weighted values as raw SUMS. */
        if (wsum <= 0.00001 && getenv("MOLE_BRANCH_AUDIT"))
            fprintf(stderr, "BRANCH weighted-raw-sums wsum=%g\n", wsum);
        if (wsum > 0.00001) {
            w.hydropathy /= wsum; w.hydrophobicity /= wsum; w.polarity /= wsum;
            w.logp /= wsum; w.logd /= wsum; w.logs /= wsum; wm /= wsum;
        }
        w.mutability = (int)wm;
    }
    w.charge = out->props.charge; w.ionizable = out->props.ionizable;
    w.npos = out->props.npos; w.nneg = out->props.nneg;
    out->wprops = w;

    out->nl = nl; out->layer = L;
    free(S); free(sd); free(sr); free(sf); free(sb); free(allres); free(allbb);
    return 0;
fail:
    free(S); free(sd); free(sr); free(sf); free(sb); free(L);
    free(allres); free(allbb); free(out->flow);
    memset(out, 0, sizeof *out);
    return -1;
}

void mole_lining_free(mole_lining *L)
{
    free(L->layer); free(L->flow);
    memset(L, 0, sizeof *L);
}

/* ---- cavity residues and properties ------------------------------------ */

/* CalculateHydrophibilicyPolarityHydropathy's RESIDUES overload
   (PhysicoChemicalProperties.cs:154). Differs from the layers one: EVERY
   residue contributes the backbone constants as well as its own side-chain
   values, and the second loop has no table lookup at all, so an unknown residue
   still adds to the denominator. */
static void hydro_block_residues(const int *res, int nres,
                                 const mole_residues *R, mole_props *P)
{
    int i, count = 0;
    P->hydropathy = P->hydrophobicity = P->polarity = 0.0;
    P->logp = P->logd = P->logs = 0.0;
    for (i = 0; i < nres; i++) {
        const res_info *info = lookup(R->name[res[i]]);
        if (!info) continue;
        count++;
        P->hydropathy += info->hydropathy;
        P->hydrophobicity += info->hydrophobicity;
        P->polarity += info->polarity;
        P->logp += info->logp; P->logd += info->logd; P->logs += info->logs;
    }
    for (i = 0; i < nres; i++) {
        count++;
        P->polarity += BB_POLARITY;
        P->hydrophobicity += BB_HYDROPHOBICITY;
        P->hydropathy += BB_HYDROPATHY;
        P->logp += BB_LOGP; P->logd += BB_LOGD; P->logs += BB_LOGS;
    }
    if (count == 0) {
        P->hydropathy = P->hydrophobicity = P->polarity = 0.0;
        P->logp = P->logd = P->logs = 0.0;
    } else {
        P->hydropathy /= count; P->hydrophobicity /= count; P->polarity /= count;
        P->logp /= count; P->logd /= count; P->logs /= count;
    }
}

void mole_cavity_properties(const int *res, int nres, const mole_residues *R,
                            mole_props *P)
{
    /* count_block skips backbone-flagged entries; a cavity residue list has no
       such split, so every entry counts. */
    static char none[1] = { 0 };
    char *flags = calloc((size_t)(nres > 0 ? nres : 1), 1);
    (void)none;
    if (!flags) { memset(P, 0, sizeof *P); return; }
    count_block(res, flags, nres, R, P);
    hydro_block_residues(res, nres, R, P);
    free(flags);
}

/* Cavity.Create's residue split. A cavity tetrahedron with fewer than four
   neighbours IN THE CAVITY contributes its unshared facets; those facets are
   boundary only when the TETRAHEDRON is (Facet.Boundary copies t.IsBoundary),
   and their three atoms' residues form BoundaryResidues. InnerResidues is every
   other residue touched by the cavity. Both sorted by chain then number. */
int mole_cavity_residues(const mole_complex *M, int comp, const int *pres,
                         const mole_residues *R, int *bnd, int *inner,
                         int *out_nb, int *out_ni)
{
    int t, k, j, nb = 0, ni = 0;
    char *isb = calloc((size_t)R->n, 1), *any = calloc((size_t)R->n, 1);
    if (!isb || !any) { free(isb); free(any); return -1; }
    for (t = 0; t < M->nt; t++) {
        int deg = 0;
        if (!M->alive[t] || M->comp[t] != comp) continue;
        for (k = 0; k < 4; k++) {
            int n = M->tn[4*t+k];
            if (n >= 0 && M->alive[n] && M->comp[n] == comp) deg++;
        }
        for (j = 0; j < 4; j++) {
            int r = pres[M->tv[4*t+j]];
            if (r >= 0) any[r] = 1;
        }
        if (deg >= 4 || !M->boundary[t]) continue;
        for (k = 0; k < 4; k++) {
            int n = M->tn[4*t+k];
            if (n >= 0 && M->alive[n] && M->comp[n] == comp) continue;
            for (j = 0; j < 4; j++) {          /* the facet omits vertex k */
                int r;
                if (j == k) continue;
                r = pres[M->tv[4*t+j]];
                if (r >= 0) isb[r] = 1;
            }
        }
    }
    for (j = 0; j < R->n; j++) {
        if (isb[j]) bnd[nb++] = j;
        else if (any[j]) inner[ni++] = j;
    }
    free(isb); free(any);
    /* Counts go out through the caller's pointers. The old `(nb<<16)|ni`
       assumed each set stays under 32768, which nothing enforces - R->n is
       bounded only by MOLE_MAXRES (60000). Past 32768 boundary residues the
       packed value reaches the sign bit, and the caller's `if (packed < 0)`
       error test then makes the whole cavity VANISH from the output. */
    if (out_nb) *out_nb = nb;
    if (out_ni) *out_ni = ni;
    return 0;
}

/* ---- Tunnel.FindHetResidues (Tunnel.cs:638) ---------------------------- */

/* Tunnel.FindHetResidues - the HET residues reported in <HetResidues>. Union of
   the atoms within 1.2 * Radius of a density-1 profile point and the vertices of
   every tetrahedron on the path, filtered to IsHetAtom and not water, then
   ordered by chain identifier and residue number.

   The search uses Structure.InvariantKdAtomTree(), keyed on InvariantPosition -
   coordinates BEFORE the jitter - while the profile points come from the
   jittered pipeline; A->invpiv holds the former. The tree covers all of
   Structure.Atoms with the water filter applied to the result, so searching
   pivots only is equivalent. Density is 1, not the default 8. */
int mole_het_residues(const mole_tunnel_profile *pr, const mole_atoms *A,
                      const double *piv, int np,
                      const mole_complex *M, const int *path, int npath,
                      const int *pres, const mole_residues *R, int *out, int cap)
{
    int n = (int)pr->length, i, k, nout = 0;
    double dt;
    char *seen = calloc((size_t)(R->n ? R->n : 1), 1);
    if (!seen) return 0;
    if (n < 1) n = 1;
    dt = 1.0 / n;
    (void)piv;

    for (i = 0; i <= n; i++) {
        double t = dt * i;
        double cx = mole_spline_eval(&pr->sx, t);
        double cy = mole_spline_eval(&pr->sy, t);
        double cz = mole_spline_eval(&pr->sz, t);
        double rr = 1.2 * mole_spline_eval(&pr->sr, t);
        double r2 = rr * rr;
        int a;
        /* MOLE_HET_PROBE: the evidence for a claim that looks like a bug.
           Radius is the distance to the atom SURFACE and the tree is keyed on
           atom CENTRES, so a centre sits at Radius + vdW and 1.2 * Radius
           reaches one only when Radius >= 5 * vdW, about 8 A. Real tunnels are
           1-3 A, so this half of the union never fires - in MOLE either, since
           the formula and the tree are both theirs. Do NOT "fix" it by adding
           vdW: that would be a deviation, and the path half below is what
           actually produces <HetResidues>. */
        if (getenv("MOLE_HET_PROBE")) {
            double best = 1e300; int q;
            for (q = 0; q < np; q++) {
                double ax = A->invpiv[3*q]-cx, ay = A->invpiv[3*q+1]-cy, az = A->invpiv[3*q+2]-cz;
                double d = sqrt(ax*ax+ay*ay+az*az);
                if (d < best) best = d;
            }
            fprintf(stderr, "HETPROBE t=%.3f R=%.3f 1.2R=%.3f nearestAtomCentre=%.3f %s\n",
                    t, rr/1.2, rr, best, best <= rr ? "REACHES" : "-");
        }
        if (rr <= 0) continue;
        for (a = 0; a < np; a++) {
            double dx = A->invpiv[3*a]   - cx;
            double dy = A->invpiv[3*a+1] - cy;
            double dz = A->invpiv[3*a+2] - cz;
            int ri;
            if (dx*dx + dy*dy + dz*dz > r2) continue;
            ri = pres[a];
            if (ri < 0 || seen[ri]) continue;
            if (!mole_is_het_resname(R->name[ri])) continue;
            seen[ri] = 1;
            if (nout < cap) out[nout++] = ri;
        }
    }
    /* ...Concat(this.Path.SelectMany(t => t.Vertices...)) - the path's own
       tetrahedra, which reach residues no profile sphere happens to cover. */
    for (i = 0; i < npath; i++)
        for (k = 0; k < 4; k++) {
            int a = M->tv[4*path[i]+k], ri;
            if (a < 0 || a >= np) continue;
            ri = pres[a];
            if (ri < 0 || seen[ri]) continue;
            if (!mole_is_het_resname(R->name[ri])) continue;
            seen[ri] = 1;
            if (nout < cap) out[nout++] = ri;
        }
    free(seen);

    mole_sort_residues_by_chain(out, nout, R);
    return nout;
}

/* OrderBy(ChainIdentifier).ThenBy(Number), used verbatim by HetResidues above.
   LINQ's sort is stable, so whatever tie the caller's Distinct order left is
   kept - insertion sort preserves that. Exposed so the chain-ordering
   regression (mole_lining_regression_test.c) can drive it with a synthetic
   residue table instead of needing full profile/path geometry. */
void mole_sort_residues_by_chain(int *out, int nout, const mole_residues *R)
{
    int i, j;
    for (i = 1; i < nout; i++) {
        int v = out[i];
        j = i - 1;
        while (j >= 0) {
            int c = mole_chain_cmp(R->chain[out[j]], R->chain[v]);
            if (c > 0 || (c == 0 && R->seq[out[j]] > R->seq[v])) { out[j+1] = out[j]; j--; }
            else break;
        }
        out[j+1] = v;
    }
}
