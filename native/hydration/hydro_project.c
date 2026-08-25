/* hydro_project.c -- compiled accelerator for VMDHole's Hydration analysis
 * (compute_hydration in vmdhole.tcl). Ports ONLY the measured hot path:
 * per-frame water residue-COG reduction, axial projection onto the pore
 * axis, and the envelope-radius in/out test that produces "qco" (the
 * per-frame list of accepted axial water coordinates) -- optionally also
 * the KDE/hard-bin accumulation that follows it (--bin).
 *
 * WHY THIS PATH: profiled on a real 50-frame solvated CHARMM trajectory
 * (step5_assembly.hmr/sim_1.dcd, 30600 waters, ~2213 bbox-candidate waters
 * per frame). Under the plugin's SHIPPED DEFAULT settings (water_kde_bw=1.4,
 * fixed -- NOT "auto"/AMISE), this projection+envelope-test loop is ~60-85%
 * of compute_hydration's total wall time and the ONLY phase that scales with
 * both frame count and water count; the AMISE bandwidth solve
 * (_amise_optimal_bandwidth/_amise_bandwidths_parallel) is NOT invoked at
 * all under the default fixed bandwidth and was NOT the bottleneck measured
 * here -- see native/NOTES-hydration-accel.md for the full phase
 * breakdown and the commands that produced it.
 *
 * BIT-IDENTICALITY: every arithmetic expression below is transcribed
 * expression-for-expression from the Tcl source (envelope_radius and the
 * qco projection/inclusion test in compute_hydration, vmdhole.tcl), same
 * operation order, same associativity. Compiled with -ffp-contract=off (see
 * native/build.sh's hydro step) so the compiler cannot fuse multiply-add pairs
 * (dxc*ux + dyc*uy + dzc*uz, ra + f*(rb-ra), ...) into a single FMA
 * instruction, which would round once instead of twice like Tcl's
 * expr-by-expr evaluation and silently break bit-identity at the last ULP
 * (this project has been bitten by exactly this class of bug before, see
 * sos_triangle_fast.c's own ULP-divergence note near hydro_thin_spheres).
 * strtod() is used for all numeric parsing (locale-independent, correctly
 * rounded) rather than scanf's %f/%lf.
 *
 * --bin's GLOBAL accumulator (see below) matters because floating-point
 * addition is not associative: Tcl's ORIGINAL code adds every individual
 * water/window contribution straight into one running global total, frame
 * by frame in `result_frames` order -- T = (((0+a1)+a2)+b1)+b2. Computing a
 * per-frame SUBTOTAL first and merging subtotals (T' = (a1+a2)+(b1+b2)) is
 * mathematically the same value but NOT bit-identical in general (e.g.
 * a1+a2 that rounds to exactly 1.0 vs. accumulating b1,b2 first changes
 * which roundings happen). frame_bincount (the PER-FRAME table) genuinely
 * is a from-zero running total local to that frame in the original code, so
 * a per-frame subtotal IS correct there -- only the cross-frame GLOBAL
 * total needs the true running-accumulation treatment below.
 *
 * CLI / I/O convention: mirrors sos_triangle_fast.c's own --batch (see its
 * "--batch" handling): a batch file lists TAB-separated "IN\tOUT" pairs,
 * one per frame; stdin/stdout are freopen'd to each pair in turn so the
 * single-job code path (stdin/stdout) and the batch path share the same
 * per-frame parser. --bin additionally requires --bin-global PATH (batch
 * mode only) to receive the ONE cross-frame global bin table, written once
 * after every frame in the batch has been processed, in the SAME frame
 * order the batch file lists them. See NOTES-hydration-accel.md for the
 * exact per-frame text format and why single-job mode does not support
 * --bin (there is nothing to distinguish "global" from "local" with one
 * frame; use plain project-only mode, or --batch with a single-line batch
 * file, if a global total is still wanted for exactly one frame).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../xalloc.h"

/* ---------------------------------------------------------------------- */
/* growable double/int arrays                                              */
/* ---------------------------------------------------------------------- */
typedef struct { double *d; int n, cap; } DArr;
typedef struct { long   *d; int n, cap; } LArr;

static void darr_push(DArr *a, double v) {
    if (a->n >= a->cap) {
        a->cap = a->cap ? a->cap * 2 : 256;
        a->d = xa_realloc(a->d, (size_t)a->cap * sizeof(double));
        if (!a->d) { fprintf(stderr, "hydro_project: out of memory\n"); exit(1); }
    }
    a->d[a->n++] = v;
}
static void larr_push(LArr *a, long v) {
    if (a->n >= a->cap) {
        a->cap = a->cap ? a->cap * 2 : 256;
        a->d = xa_realloc(a->d, (size_t)a->cap * sizeof(long));
        if (!a->d) { fprintf(stderr, "hydro_project: out of memory\n"); exit(1); }
    }
    a->d[a->n++] = v;
}

/* ---------------------------------------------------------------------- */
/* line reading: one logical line, arbitrary length, no trailing newline   */
/* ---------------------------------------------------------------------- */
static char *read_line(FILE *f) {
    size_t cap = 4096, len = 0;
    char *buf = malloc(cap);
    if (!buf) { fprintf(stderr, "hydro_project: out of memory\n"); exit(1); }
    int c;
    int any = 0;
    while ((c = fgetc(f)) != EOF) {
        any = 1;
        if (c == '\n') break;
        if (len + 1 >= cap) { cap *= 2; buf = realloc(buf, cap); if (!buf) { fprintf(stderr,"hydro_project: OOM\n"); exit(1);} }
        buf[len++] = (char)c;
    }
    if (!any && c == EOF) { free(buf); return NULL; }
    buf[len] = 0;
    return buf;
}

static void parse_doubles_line(const char *line, int n, DArr *out) {
    const char *p = line;
    int i;
    for (i = 0; i < n; i++) {
        char *end;
        double v = strtod(p, &end);
        if (end == p) { fprintf(stderr, "hydro_project: expected %d doubles, got %d\n", n, i); exit(1); }
        darr_push(out, v);
        p = end;
    }
}
static void parse_longs_line(const char *line, int n, LArr *out) {
    const char *p = line;
    int i;
    for (i = 0; i < n; i++) {
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) { fprintf(stderr, "hydro_project: expected %d ints, got %d\n", n, i); exit(1); }
        larr_push(out, v);
        p = end;
    }
}

/* ---------------------------------------------------------------------- */
/* envelope_radius: EXACT port of ::VMDHole::envelope_radius (vmdhole.tcl). */
/* Linear scan, first bracketing interval wins (env may contain duplicate  */
/* co values -- do NOT replace with a binary search, see NOTES).           */
/* ---------------------------------------------------------------------- */
static double envelope_radius(const double *ec, const double *er, int n, double co) {
    int i;
    if (n == 0) return 0.0;
    if (co <= ec[0])     return er[0];
    if (co >= ec[n - 1]) return er[n - 1];
    for (i = 0; i < n - 1; i++) {
        double ca = ec[i], cb = ec[i + 1];
        if (co >= ca && co <= cb) {
            double ra = er[i], rb = er[i + 1];
            double f = (cb > ca) ? (co - ca) / (cb - ca) : 0.0;
            return ra + f * (rb - ra);
        }
    }
    return er[n - 1];
}

/* ---------------------------------------------------------------------- */
/* Residue centre-of-geometry reduction: EXACT port of the reduction block */
/* in compute_hydration (vmdhole.tcl) -- one output position per DISTINCT  */
/* resid, in FIRST-ENCOUNTER order (matches Tcl's `dict keys $_wcn`,       */
/* insertion-ordered). A simple open-addressing hash keyed on resid gives  */
/* O(1) amortised lookup while preserving that order via a parallel        */
/* first-seen index array.                                                 */
/* ---------------------------------------------------------------------- */
typedef struct {
    long  resid;
    double sx, sy, sz;
    long  cnt;
    int   used;
} CogSlot;

static int cog_reduce(const long *resid, const double *x, const double *y, const double *z,
                       int n, DArr *ox, DArr *oy, DArr *oz) {
    if (n == 0) return 0;
    size_t cap = 16;
    while (cap < (size_t)n * 4) cap *= 2;
    CogSlot *tab = calloc(cap, sizeof(CogSlot));
    if (!tab) { fprintf(stderr, "hydro_project: out of memory\n"); exit(1); }
    int *order = malloc((size_t)n * sizeof(int));
    if (!order) { fprintf(stderr, "hydro_project: out of memory\n"); exit(1); }
    int norder = 0;
    int i;
    for (i = 0; i < n; i++) {
        long rid = resid[i];
        size_t h = ((size_t)(rid) * 2654435761u) & (cap - 1);
        while (tab[h].used && tab[h].resid != rid) h = (h + 1) & (cap - 1);
        if (!tab[h].used) {
            tab[h].used = 1;
            tab[h].resid = rid;
            tab[h].sx = tab[h].sy = tab[h].sz = 0.0;
            tab[h].cnt = 0;
            order[norder++] = (int)h;
        }
        tab[h].sx += x[i];
        tab[h].sy += y[i];
        tab[h].sz += z[i];
        tab[h].cnt += 1;
    }
    for (i = 0; i < norder; i++) {
        CogSlot *s = &tab[order[i]];
        darr_push(ox, s->sx / (double)s->cnt);
        darr_push(oy, s->sy / (double)s->cnt);
        darr_push(oz, s->sz / (double)s->cnt);
    }
    free(tab);
    free(order);
    return norder;
}

/* ---------------------------------------------------------------------- */
/* Per-frame job: parsed input + derived qco.                              */
/* ---------------------------------------------------------------------- */
typedef struct {
    char  *out_path;
    DArr   qco;
    double dz, bw;
    int    use_kde;
    int    have_bin_params;
    long   n_raw_waters; /* pre-COG atom count (for nfwater semantics, see NOTES) */
} FrameJob;

static void fail(const char *msg) { fprintf(stderr, "hydro_project: %s\n", msg); exit(1); }

/* Parse one job from the CURRENT stdin and compute its qco. Does not write
   any output (that happens later, see write_frame_output / global pass). */
static void parse_and_project(FrameJob *job) {
    char *line;
    double mx=0,my=0,mz=0,ux=0,uy=0,uz=0, cmin=0,cmax=0, dcap=0;
    DArr env_co = {0}, env_r = {0};
    LArr w_resid = {0};
    DArr w_x = {0}, w_y = {0}, w_z = {0};
    int n_env = 0, n_w = 0;

    line = read_line(stdin);
    if (!line || strncmp(line, "AXIS", 4) != 0) fail("expected AXIS line");
    if (sscanf(line + 4, "%lf %lf %lf %lf %lf %lf", &mx,&my,&mz,&ux,&uy,&uz) != 6)
        fail("malformed AXIS line");
    free(line);

    line = read_line(stdin);
    if (!line || strncmp(line, "RANGE", 5) != 0) fail("expected RANGE line");
    if (sscanf(line + 5, "%lf %lf", &cmin, &cmax) != 2) fail("malformed RANGE line");
    free(line);

    line = read_line(stdin);
    if (!line || strncmp(line, "DCAP", 4) != 0) fail("expected DCAP line");
    if (sscanf(line + 4, "%lf", &dcap) != 1) fail("malformed DCAP line");
    free(line);

    line = read_line(stdin);
    if (!line || strncmp(line, "ENV", 3) != 0) fail("expected ENV line");
    if (sscanf(line + 3, "%d", &n_env) != 1) fail("malformed ENV line");
    free(line);
    if (n_env > 0) {
        line = read_line(stdin); if (!line) fail("missing ENV coords line");
        parse_doubles_line(line, n_env, &env_co); free(line);
        line = read_line(stdin); if (!line) fail("missing ENV radii line");
        parse_doubles_line(line, n_env, &env_r); free(line);
    }

    line = read_line(stdin);
    if (!line || strncmp(line, "WATERS", 6) != 0) fail("expected WATERS line");
    if (sscanf(line + 6, "%d", &n_w) != 1) fail("malformed WATERS line");
    free(line);
    if (n_w > 0) {
        line = read_line(stdin); if (!line) fail("missing WATERS resid line");
        parse_longs_line(line, n_w, &w_resid); free(line);
        line = read_line(stdin); if (!line) fail("missing WATERS x line");
        parse_doubles_line(line, n_w, &w_x); free(line);
        line = read_line(stdin); if (!line) fail("missing WATERS y line");
        parse_doubles_line(line, n_w, &w_y); free(line);
        line = read_line(stdin); if (!line) fail("missing WATERS z line");
        parse_doubles_line(line, n_w, &w_z); free(line);
    }
    job->n_raw_waters = n_w;

    job->have_bin_params = 0;
    job->dz = 0; job->bw = 0; job->use_kde = 0;
    {
        long save = ftell(stdin);
        line = read_line(stdin);
        if (line && strncmp(line, "BIN", 3) == 0) {
            if (sscanf(line + 3, "%lf %lf %d", &job->dz, &job->bw, &job->use_kde) != 3) fail("malformed BIN line");
            job->have_bin_params = 1;
            free(line);
        } else {
            if (line) free(line);
            fseek(stdin, save, SEEK_SET);
        }
    }

    /* residue COG reduction (always performed - matches Tcl's near-universal
       branch; see NOTES-hydration-accel.md for the one unreachable-in-
       practice Tcl branch this does not replicate). */
    DArr rx = {0}, ry = {0}, rz = {0};
    int n_res = cog_reduce(w_resid.d, w_x.d, w_y.d, w_z.d, n_w, &rx, &ry, &rz);

    job->qco.d = NULL; job->qco.n = 0; job->qco.cap = 0;
    int i;
    for (i = 0; i < n_res; i++) {
        double wx = rx.d[i], wy = ry.d[i], wz = rz.d[i];
        double dxc = wx - mx, dyc = wy - my, dzc = wz - mz;
        double co = dxc*ux + dyc*uy + dzc*uz;
        if (co < cmin || co > cmax) continue;
        double rr = envelope_radius(env_co.d, env_r.d, n_env, co);
        if (rr <= 0) continue;
        double rr_eff = (dcap > 0 && rr > dcap) ? dcap : rr;
        double perp2 = dxc*dxc + dyc*dyc + dzc*dzc - co*co;
        if (perp2 < 0) perp2 = 0;
        if (perp2 <= rr_eff * rr_eff) darr_push(&job->qco, co);
    }

    free(env_co.d); free(env_r.d);
    free(w_resid.d); free(w_x.d); free(w_y.d); free(w_z.d);
    free(rx.d); free(ry.d); free(rz.d);
}

static void frame_bin_range(const FrameJob *job, long *lo, long *hi) {
    int i;
    int have = 0;
    for (i = 0; i < job->qco.n; i++) {
        double co = job->qco.d[i];
        long bi_lo, bi_hi;
        if (job->use_kde) {
            bi_lo = (long)floor((co - 3.0*job->bw) / job->dz);
            bi_hi = (long)floor((co + 3.0*job->bw) / job->dz);
        } else {
            bi_lo = bi_hi = (long)floor(co / job->dz);
        }
        if (!have) { *lo = bi_lo; *hi = bi_hi; have = 1; }
        else { if (bi_lo < *lo) *lo = bi_lo; if (bi_hi > *hi) *hi = bi_hi; }
    }
    if (!have) { *lo = 0; *hi = -1; }
}

static void write_qco_and_bins(FILE *out, const FrameJob *job, const double *local_bins, long lo, long hi) {
    int i;
    fprintf(out, "QCO %d\n", job->qco.n);
    for (i = 0; i < job->qco.n; i++) fprintf(out, "%s%.17g", i ? " " : "", job->qco.d[i]);
    fprintf(out, "\n");
    if (job->have_bin_params) {
        long nbins = hi - lo + 1; if (nbins < 0) nbins = 0;
        fprintf(out, "BINS %ld %ld\n", lo, hi);
        for (i = 0; i < nbins; i++) fprintf(out, "%s%.17g", i ? " " : "", local_bins[i]);
        fprintf(out, "\n");
    }
}

/* ---------------------------------------------------------------------- */
static int g_bin_mode = 0;      /* --bin */
static const char *g_bin_global_path = NULL; /* --bin-global PATH (batch only) */

int main(int argc, char **argv) {
    int ai;
    const char *batch_file = NULL;

    for (ai = 1; ai < argc; ai++) {
        if (strcmp(argv[ai], "--hole-features") == 0) {
            printf("hole_features: hydroproject hydrobin\n");
            return 0;
        } else if (strcmp(argv[ai], "--bin") == 0) {
            g_bin_mode = 1;
        } else if (strcmp(argv[ai], "--bin-global") == 0) {
            if (ai + 1 >= argc) { fprintf(stderr, "hydro_project: --bin-global requires a FILE argument\n"); return 1; }
            g_bin_global_path = argv[++ai];
        } else if (strcmp(argv[ai], "--batch") == 0) {
            if (ai + 1 >= argc) { fprintf(stderr, "hydro_project: --batch requires a FILE argument\n"); return 1; }
            batch_file = argv[++ai];
        } else {
            fprintf(stderr, "hydro_project: unrecognized argument: %s\n", argv[ai]);
            return 1;
        }
    }

    if (!batch_file) {
        /* single-job mode: stdin -> stdout. No global accumulator (nothing
           to distinguish it from the local one with a single frame). */
        if (g_bin_mode && g_bin_global_path) fail("--bin-global requires --batch");
        FrameJob job;
        memset(&job, 0, sizeof(job));
        parse_and_project(&job);
        if (g_bin_mode && !job.have_bin_params) fail("--bin given but job has no BIN line");
        double *local_bins = NULL;
        long lo = 0, hi = -1;
        if (g_bin_mode) {
            frame_bin_range(&job, &lo, &hi);
            long nbins = hi - lo + 1; if (nbins < 0) nbins = 0;
            local_bins = nbins > 0 ? xa_calloc((size_t)nbins, sizeof(double)) : NULL;
            int i;
            for (i = 0; i < job.qco.n; i++) {
                double co = job.qco.d[i];
                if (job.use_kde) {
                    double norm = job.dz / (job.bw * 2.5066282746310002);
                    long blo = (long)floor((co - 3.0*job.bw) / job.dz);
                    long bhi = (long)floor((co + 3.0*job.bw) / job.dz);
                    long bi;
                    for (bi = blo; bi <= bhi; bi++) {
                        double dzz = ((double)bi + 0.5) * job.dz - co;
                        double wkde = norm * exp(-dzz*dzz / (2.0*job.bw*job.bw));
                        local_bins[bi - lo] += wkde;
                    }
                } else {
                    long bi = (long)floor(co / job.dz);
                    local_bins[bi - lo] += 1.0;
                }
            }
        }
        write_qco_and_bins(stdout, &job, local_bins, lo, hi);
        free(local_bins);
        free(job.qco.d);
        return 0;
    }

    /* --batch mode ------------------------------------------------------ */
    if (g_bin_mode && !g_bin_global_path) fail("--bin --batch requires --bin-global PATH");

    FILE *bf = fopen(batch_file, "r");
    if (!bf) { fprintf(stderr, "hydro_project: cannot open batch file: %s\n", batch_file); return 1; }

    /* Pass 1: parse every job and compute its qco (kept in memory - total
       size is O(total qualifying waters across the whole trajectory), not
       O(all candidate waters), so this is small even at 10k+ frames). */
    FrameJob *jobs = NULL; int njobs = 0, cap_jobs = 0;
    char line[8192];
    while (fgets(line, sizeof(line), bf)) {
        char in_path[4096], out_path[4096];
        char *tok;
        line[strcspn(line, "\n")] = 0;
        if (line[0] == '#' || line[0] == '\0') continue;
        tok = strtok(line, "\t");
        if (!tok) continue;
        strncpy(in_path, tok, sizeof(in_path)-1); in_path[sizeof(in_path)-1]=0;
        tok = strtok(NULL, "\t");
        if (!tok) continue;
        strncpy(out_path, tok, sizeof(out_path)-1); out_path[sizeof(out_path)-1]=0;
        if (freopen(in_path, "r", stdin) == NULL) {
            fprintf(stderr, "hydro_project: cannot open input: %s\n", in_path);
            continue;
        }
        if (njobs >= cap_jobs) {
            cap_jobs = cap_jobs ? cap_jobs * 2 : 64;
            jobs = realloc(jobs, (size_t)cap_jobs * sizeof(FrameJob));
            if (!jobs) fail("out of memory (jobs array)");
        }
        FrameJob *job = &jobs[njobs++];
        memset(job, 0, sizeof(*job));
        job->out_path = xa_strdup(out_path);
        parse_and_project(job);
        if (g_bin_mode && !job->have_bin_params) fail("--bin given but a job has no BIN line");
    }
    fclose(bf);

    /* Pass 2: write each frame's own (local) QCO+BINS, and -- in the SAME
       frame order, contribution by contribution -- accumulate the ONE
       cross-frame GLOBAL bin table exactly as Tcl's original nested loop
       would (see header comment: this is the part that must NOT be done
       as a per-frame subtotal merge). */
    long glo = 0, ghi = -1;
    int have_global = 0;
    if (g_bin_mode) {
        int f;
        for (f = 0; f < njobs; f++) {
            long lo, hi;
            frame_bin_range(&jobs[f], &lo, &hi);
            if (jobs[f].qco.n == 0) continue;
            if (!have_global) { glo = lo; ghi = hi; have_global = 1; }
            else { if (lo < glo) glo = lo; if (hi > ghi) ghi = hi; }
        }
    }
    double *global_bins = NULL;
    if (g_bin_mode && have_global) {
        long gn = ghi - glo + 1;
        global_bins = calloc((size_t)gn, sizeof(double));
        if (!global_bins) fail("out of memory (global bins)");
    }

    int f;
    for (f = 0; f < njobs; f++) {
        FrameJob *job = &jobs[f];
        double *local_bins = NULL;
        long lo = 0, hi = -1;
        if (g_bin_mode) {
            frame_bin_range(job, &lo, &hi);
            long nbins = hi - lo + 1; if (nbins < 0) nbins = 0;
            local_bins = nbins > 0 ? xa_calloc((size_t)nbins, sizeof(double)) : NULL;
            int i;
            for (i = 0; i < job->qco.n; i++) {
                double co = job->qco.d[i];
                if (job->use_kde) {
                    double norm = job->dz / (job->bw * 2.5066282746310002);
                    long blo = (long)floor((co - 3.0*job->bw) / job->dz);
                    long bhi = (long)floor((co + 3.0*job->bw) / job->dz);
                    long bi;
                    for (bi = blo; bi <= bhi; bi++) {
                        double dzz = ((double)bi + 0.5) * job->dz - co;
                        double wkde = norm * exp(-dzz*dzz / (2.0*job->bw*job->bw));
                        local_bins[bi - lo] += wkde;
                        if (global_bins) global_bins[bi - glo] += wkde;
                    }
                } else {
                    long bi = (long)floor(co / job->dz);
                    local_bins[bi - lo] += 1.0;
                    if (global_bins) global_bins[bi - glo] += 1.0;
                }
            }
        }
        FILE *out = fopen(job->out_path, "w");
        if (!out) { fprintf(stderr, "hydro_project: cannot open output: %s\n", job->out_path); }
        else {
            write_qco_and_bins(out, job, local_bins, lo, hi);
            fclose(out);
        }
        free(local_bins);
    }

    if (g_bin_mode) {
        FILE *gf = fopen(g_bin_global_path, "w");
        if (!gf) fail("cannot open --bin-global output path");
        if (have_global) {
            long gn = ghi - glo + 1;
            int i;
            fprintf(gf, "BINS %ld %ld\n", glo, ghi);
            for (i = 0; i < gn; i++) fprintf(gf, "%s%.17g", i ? " " : "", global_bins[i]);
            fprintf(gf, "\n");
        } else {
            fprintf(gf, "BINS 0 -1\n\n");
        }
        fclose(gf);
    }

    for (f = 0; f < njobs; f++) { free(jobs[f].out_path); free(jobs[f].qco.d); }
    free(jobs);
    free(global_bins);
    return 0;
}
