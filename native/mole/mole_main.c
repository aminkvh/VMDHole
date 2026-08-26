/* --tunnel-mole: MOLE 2's tunnel algorithm, as an entry point.
 *
 * A third engine alongside --tunnel-search (lattice) and --tunnel-voronoi
 * (ours). This one is a port of MOLE 2 and is validated against the reference
 * implementation rather than against a principle: on 1tqn, 1BL8, 2ACE and 1AKD
 * every value of every profile row matches MOLE's own output - 2252 rows across
 * 17 tunnels.
 *
 * Output is the same "T id bottleneck length cost throughput" / "P id x y z r"
 * shape the other two engines emit, so callers do not need a third parser.
 *
 * A SEPARATE binary rather than a subcommand of sos_triangle_fast, and not for
 * convenience: this engine needs the predicates compiled at VP_SCALE 1e5 so
 * MOLE's +-0.00005 A jitter survives quantisation, while sos_triangle_fast and
 * --tunnel-voronoi need the default 1e3 and have byte-identity tests pinned to
 * it. One executable cannot hold both, and changing the shared default to suit
 * this engine would silently move every existing result.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mole_lining.h"
#include "mole_dh.h"
#include "../xalloc.h"

static mole_atoms MA;
static mole_residues MR;
static double mpiv[3*MOLE_MAXATOM], mrad[MOLE_MAXATOM], mbfac[MOLE_MAXATOM];
static int mfree[MOLE_MAXATOM], mbb[MOLE_MAXATOM], mres[MOLE_MAXATOM];

/* `path` is retained, not just its length: FindHetResidues unions the profile
   sphere search with the tunnel PATH's own tetrahedron vertices. */
typedef struct { double length, bottleneck; int plen, cav; int *path, npath;
                 mole_tunnel_profile prof; } mtun;

/* TunnelCollection.TunnelComparer: (Cavity.Id, Length) lexicographically, NOT
   length alone. Cavities are visited in ascending index here and MOLE numbers
   them the same way, so ordering on our component index reproduces theirs. The
   difference only shows when one run yields tunnels from more than one cavity -
   a user origin that snaps into two, or auto origins - which is why every
   single-cavity reference agreed under a plain length sort. */
static int by_cav_len(const void *a, const void *b)
{
    const mtun *x = a, *y = b;
    double d;
    if (x->cav != y->cav) return x->cav < y->cav ? -1 : 1;
    d = x->length - y->length;
    return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

/* Widths at %.5f and properties at %.4f: MOLE's own XML prints the widths to 5
   decimals, so a comparison against it is limited by our digits otherwise. */
static void write_props(FILE *out, const char *tag, int id, const mole_props *P)
{
    fprintf(out, "%s %d %d %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %d\n",
            tag, id, P->charge, P->ionizable, P->npos, P->nneg,
            P->hydropathy, P->hydrophobicity, P->polarity,
            P->logp, P->logd, P->logs, P->mutability);
}

/* <HetResidues>. Separate from the lining: different search radius
   (1.2 * Radius, not the nearest five), un-jittered coordinates, no backbone
   split. */
static void write_het(FILE *out, int id, const mole_tunnel_profile *pr,
                      const mole_complex *M, const int *path, int npath, int np)
{
    int het[MOLE_MAXRES], n, j;
    n = mole_het_residues(pr, &MA, mpiv, np, M, path, npath, mres, &MR,
                          het, MOLE_MAXRES);
    fprintf(out, "H %d %d", id, n);
    for (j = 0; j < n; j++)
        fprintf(out, " %s:%d:%s", MR.name[het[j]], MR.seq[het[j]], MR.chain[het[j]]);
    fputc('\n', out);
}

static void write_lining(FILE *out, int id, const mole_tunnel_profile *pr, int np)
{
    mole_lining L;
    int j, k;
    if (mole_lining_build(pr, mpiv, mrad, np, mres, &MR, mbb, &L) != 0) return;
    for (j = 0; j < L.nl; j++) {
        const mole_layer *y = &L.layer[j];
        fprintf(out, "L %d %d %.5f %.5f %.5f %.5f %.5f %d %d %d %d %d"
                     " %.4f %.4f %.4f %.4f %.4f %.4f %d %d",
                id, j + 1, y->start, y->end, y->radius, y->freeradius, y->bradius,
                y->localmin, y->props.charge, y->props.ionizable,
                y->props.npos, y->props.nneg, y->props.hydropathy,
                y->props.hydrophobicity, y->props.polarity,
                y->props.logp, y->props.logd, y->props.logs,
                y->props.mutability, y->nres);
        for (k = 0; k < y->nres; k++)
            fprintf(out, " %s:%d:%s:%d:%d", MR.name[y->res[k]], MR.seq[y->res[k]],
                    MR.chain[y->res[k]], (int)y->isbb[k], y->flow[k]);
        fputc('\n', out);
    }
    fprintf(out, "F %d %d", id, L.nflow);
    for (j = 0; j < L.nflow; j++)
        fprintf(out, " %s:%d:%s:%d", MR.name[L.flow[j].res], MR.seq[L.flow[j].res],
                MR.chain[L.flow[j].res], (int)L.flow[j].isbb);
    fputc('\n', out);
    write_props(out, "Y", id, &L.props);
    write_props(out, "W", id, &L.wprops);
    mole_lining_free(&L);
}

/* One cavity's residue split and properties, as MOLE's cavities.xml reports
   them. Volume/Depth/DepthLength were already computed; the residues and the
   physico-chemical block are what this adds. */
static void write_cavity(FILE *out, int id, const mole_complex *M, int comp,
                         const mole_cavity *cv, const int *pres)
{
    int *bnd = xa_malloc((size_t)(MR.n ? MR.n : 1) * sizeof(int));
    int *inn = malloc((size_t)(MR.n ? MR.n : 1) * sizeof(int));
    int packed, nb, ni, j;
    mole_props P;
    if (!bnd || !inn) { free(bnd); free(inn); return; }
    packed = mole_cavity_residues(M, comp, pres, &MR, bnd, inn, &nb, &ni);
    if (packed < 0) { free(bnd); free(inn); return; }
    /* Type: MOLE calls a cavity with no boundary residues a Void. */
    fprintf(out, "V %d %s %.3f %d %.17g %d %d\n", id, nb ? "Cavity" : "Void",
            cv->volume, cv->depth, cv->depth_length, nb, ni);
    for (j = 0; j < 2; j++) {
        const int *set = j ? inn : bnd;
        int n = j ? ni : nb, q;
        mole_cavity_properties(set, n, &MR, &P);
        fprintf(out, "%s %d %d %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f %d %d",
                j ? "VI" : "VB", id, P.charge, P.ionizable, P.npos, P.nneg,
                P.hydropathy, P.hydrophobicity, P.polarity,
                P.logp, P.logd, P.logs, P.mutability, n);
        for (q = 0; q < n; q++)
            fprintf(out, " %s:%d:%s", MR.name[set[q]], MR.seq[set[q]],
                    MR.chain[set[q]]);
        fputc('\n', out);
    }
    free(bnd); free(inn);
}

/* GetTunnels' SurfaceCavity source, for ONE origin.
 *
 * `sources` is the origin's cavity concatenated with the SurfaceCavity
 * (TunnelComputation.cs:89) and TunnelOriginCollection.cs:106 makes one origin
 * per cavity, so the surface is a source for EVERY origin and its Dijkstra runs
 * from that origin's own tetrahedron - measured with MOLE_EXIT_DEBUG on the
 * oracle, where the printed origin differs per source block. A cavity
 * tetrahedron is always a surface member, having survived the peel and the
 * interior removal, so osrc needs no snapping.
 *
 * The SurfaceCavity has no openings of its own (UpdateOpenings returns
 * immediately for it), so it yields tunnels only through a user exit, and its
 * Cavity.Id is 0. It runs even when the origin's own cavity had NO openings -
 * MOLE prints exactly that on the 1tqn surface-exit case.
 */
static void surface_tunnel(mole_complex *Mp, const mole_params *Pp, int osrc,
                           const double *exitp, double origin_radius,
                           double bottleneck, double bottle_tol,
                           mtun **resp, int *nresp, int *ncapp, int np)
{
    mole_complex M = *Mp;
    mole_params P = *Pp;
    mtun *res = *resp;
    int nres = *nresp, ncap = *ncapp;
    {
        int pivot = mole_cavity_opening(&M, -1, M.surface, exitp, origin_radius);
        int k;
        if (getenv("MOLE_EXIT_DEBUG"))
            fprintf(stderr, "SURFACE-EXIT pivot=%d osrc=%d%s\n", pivot, osrc,
                    pivot >= 0 ? "" : "  -> skipped");
        if (pivot >= 0) {
            double *dist = xa_malloc((size_t)M.nt*sizeof(double));
            int *prev = xa_malloc((size_t)M.nt*sizeof(int));
            if (dist && prev) {
                mole_dijkstra_mask(&M, -1, osrc, dist, prev, P.weight, M.surface);
                if (getenv("MOLE_EXIT_DEBUG"))
                    fprintf(stderr, "  dijkstra dist[pivot]=%g %s\n", dist[pivot],
                            dist[pivot] < 1e299 ? "reachable" : "UNREACHABLE");
                if (dist[pivot] < 1e299) {
                    int len = 0, *p, j, ncp, *cp;
                    mole_tunnel_profile pr;
                    for (k = pivot; k >= 0; k = prev[k]) len++;
                    p = malloc((size_t)len*sizeof(int));
                    if (p) {
                            for (k = pivot, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
                        /* Tunnel.Create:680, path.TakeWhile(p => !p.IsBoundary).
                           The same truncation the cavity loop above does. It
                           matters far more here: on the surface graph the peel
                           flagged almost everything boundary, so the route is
                           cut back to the part still inside the protein. */
                        for (j = 0; j < len; j++) if (M.boundary[p[j]]) break;
                        len = j;
                        cp = len > 0 ? xa_malloc((size_t)len*sizeof(int)) : NULL;
                        ncp = cp ? mole_control_path(&M, p, len, mpiv, mrad, np,
                                                     P.interior_threshold, cp) : 0;
                        if (getenv("MOLE_EXIT_DEBUG"))
                            fprintf(stderr, "  path len=%d controlPath=%d\n", len, ncp);
                        if (ncp >= 5 && mole_profile(&M, cp, ncp, mpiv, mrad, np, &pr) == 0) {
                            if (getenv("MOLE_EXIT_DEBUG"))
                                fprintf(stderr, "  profile len=%.3f bottleneckOK=%d\n",
                                        pr.length,
                                        mole_filter_bottleneck(&pr, bottleneck, bottle_tol, 8.0));
                            if (pr.length >= P.min_tunnel_length
                                && mole_filter_bottleneck(&pr, bottleneck, bottle_tol, 8.0)) {
                                if (nres == ncap) { ncap = ncap ? ncap*2 : 32;
                                    res = xa_realloc(res, (size_t)ncap*sizeof(*res)); }
                                res[nres].length = pr.length; res[nres].plen = len;
                                res[nres].cav = 0; res[nres].prof = pr;
                                res[nres].npath = len;
                                res[nres].path = xa_malloc((size_t)(len ? len : 1) * sizeof(int));
                                if (res[nres].path) memcpy(res[nres].path, p, (size_t)len * sizeof(int));
                                else res[nres].npath = 0;
                                nres++;
                            } else mole_profile_free(&pr);
                        }
                        free(cp); free(p);
                    }
                }
            }
            free(dist); free(prev);
        }
    }
    *resp = res; *nresp = nres; *ncapp = ncap;
}

int main(int argc, char **argv)
{
    /* zeroed: with the DH triangulation the dt_mesh is never built, and the
       cleanup at the end frees it unconditionally. */
    dt_mesh m = {0};
    mole_complex M;
    /* zeroed: a field added to mole_params later defaults to MOLE's. */
    mole_params P = {0};
    mole_cavity *cav = NULL;
    int *crank = NULL;
    char *bres = NULL;      /* per-cavity boundary-residue mask, FilterBoundaryLayers */
    mtun *res = NULL;
    FILE *out;
    const char *atoms_file, *out_file;
    double ox = 0, oy = 0, oz = 0, origin_radius = 5.0;
    /* The remaining ComplexParameters, at MOLE's defaults. Taken as --key=value
       flags rather than more positional arguments: the origin block is already
       conditional, so anything appended after it shifts depending on whether a
       start point was given. */
    double surf_cover = 10.0, auto_cover = 10.0, bottleneck = 1.25;
    double bottle_tol = 0.0, max_sim = 0.9;
    int max_origins = 5, filter_boundary = 0;
    /* Custom exits: a user point that tunnels are computed TO. MOLE offers each
       to the SurfaceCavity and to every regular cavity, and with
       UseCustomExitsOnly only those user openings are used. */
    double exitp[3]; int have_exit = 0, exits_only = 0;
    /* Paths: Dijkstra between two user points on the SurfaceCavity graph.
       Complex.GetPaths - both endpoints given, no openings involved. */
    double patha[3], pathb[3]; int have_path = 0;
    int have_origin = 0, np, nc, nch, nvd, i, k, c, nres = 0, ncap = 0;

    if (argc < 3) {
        fprintf(stderr,
          "usage: %s ATOMS OUT [probe] [interior] [mindepth] [mindepthlen]\n"
          "          [mintunlen] [weight] [ox oy oz] [originradius]\n"
          "  ATOMS  \"x y z element is_water chain seq resname\" per line, FILE ORDER\n"
          "  weight 0=VoronoiScale (MOLE default) 1=LengthAndRadius 2=Length 3=Constant\n"
          "  give ox oy oz to pin the start; otherwise origins are computed\n"
          "  --cover=10 --autocover=10 --maxorigins=5 --bottleneck=1.25\n"
          "  --bottletol=0 --maxsim=0.9 --fbl=0   the rest of ComplexParameters\n"
          "  --exit=x,y,z --exitsonly=1           custom exit (CustomExits)\n"
          "  --path=x,y,z,x,y,z                   path between two points\n"
          "  --strict-interior=1                  MOLE StrictInterior; DEGENERATE\n"
          "                                       upstream - yields no cavities at all\n", argv[0]);
        return 2;
    }
    atoms_file = argv[1]; out_file = argv[2];
    /* Each positional is taken only if it is really a positional. Without the
       `--` test, `engine atoms out --cover=6` read "--cover=6" as the probe
       radius: atof gives 0.0, the run reports "0 channels, 0 voids" and exits
       0 with no diagnostic. argv[9..11] were already guarded this way below;
       argv[3..8] were not. */
#define POSN(i) ((argc > (i)) && strncmp(argv[i], "--", 2))
    P.probe_radius       = POSN(3) ? atof(argv[3]) : 3.0;
    P.interior_threshold = POSN(4) ? atof(argv[4]) : 1.25;
    P.min_depth          = POSN(5) ? atoi(argv[5]) : 8;
    P.min_depth_length   = POSN(6) ? atof(argv[6]) : 5.0;
    P.min_tunnel_length  = POSN(7) ? atof(argv[7]) : 0.0;
    P.weight             = POSN(8) ? (mole_weight_fn)atoi(argv[8]) : MOLE_W_VORONOI_SCALE;
    /* L9: an out-of-range weight silently ran MOLE's default, because
       mole_edge_cost's `default:` arm returns the same thing weight 0 selects -
       a faithful rendering of a C# switch over an enum that could not be out of
       range, but here the value comes from argv. Say so rather than pretend the
       user got what they asked for; the flag parser below already rejects an
       unknown option with exit 2. */
    if (P.weight < MOLE_W_VORONOI_SCALE || P.weight > MOLE_W_CONSTANT) {
        fprintf(stderr, "--tunnel-mole: weight %d is out of range "
                "(0=VoronoiScale 1=LengthAndRadius 2=Length 3=Constant)\n",
                (int)P.weight);
        return 2;
    }
    /* The origin is positional but the flags below are not, so a run with
       flags and no start point must not have three of them read as x y z. */
    if (argc > 11 && strncmp(argv[9], "--", 2) && strncmp(argv[10], "--", 2)
        && strncmp(argv[11], "--", 2)) {
        ox = atof(argv[9]); oy = atof(argv[10]); oz = atof(argv[11]); have_origin = 1;
    }
    if (have_origin && argc > 12 && strncmp(argv[12], "--", 2))
        origin_radius = atof(argv[12]);
    for (i = 3; i < argc; i++) {
        const char *a = argv[i], *v = strchr(a, '=');
        if (strncmp(a, "--", 2) || !v) continue;
        v++;
        if      (!strncmp(a, "--cover=", 8))       surf_cover = atof(v);
        else if (!strncmp(a, "--autocover=", 12))  auto_cover = atof(v);
        else if (!strncmp(a, "--maxorigins=", 13)) max_origins = atoi(v);
        else if (!strncmp(a, "--bottleneck=", 13)) bottleneck = atof(v);
        else if (!strncmp(a, "--bottletol=", 12))  bottle_tol = atof(v);
        else if (!strncmp(a, "--maxsim=", 9))      max_sim = atof(v);
        else if (!strncmp(a, "--fbl=", 6))         filter_boundary = atoi(v);
        else if (!strncmp(a, "--strict-interior=", 18)) P.strict_interior = atoi(v);
        else if (!strncmp(a, "--exit=", 7)) {
            if (sscanf(v, "%lf,%lf,%lf", &exitp[0], &exitp[1], &exitp[2]) == 3)
                have_exit = 1;
            else { fprintf(stderr, "--tunnel-mole: --exit needs x,y,z\n"); return 2; }
        }
        else if (!strncmp(a, "--exitsonly=", 12)) exits_only = atoi(v);
        else if (!strncmp(a, "--path=", 7)) {
            if (sscanf(v, "%lf,%lf,%lf,%lf,%lf,%lf", &patha[0], &patha[1], &patha[2],
                       &pathb[0], &pathb[1], &pathb[2]) == 6) have_path = 1;
            else { fprintf(stderr, "--tunnel-mole: --path needs x,y,z,x,y,z\n"); return 2; }
        }
        else { fprintf(stderr, "--tunnel-mole: unknown option %s\n", a); return 2; }
    }
    if (max_origins < 1) max_origins = 1;
    /* Hard upper bound: `origins` below is a fixed stack array of
       MOLE_MAX_ORIGINS. Without this an unclamped --maxorigins plus a
       cavity with that many qualifying maxima overwrites the stack. */
    if (max_origins > MOLE_MAX_ORIGINS) max_origins = MOLE_MAX_ORIGINS;

    if (mole_read_atoms(atoms_file, &MA, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "--tunnel-mole: cannot read %s\n", atoms_file); return 1;
    }
    np = mole_pivots(&MA, mpiv, mrad, NULL);
    if (np < 0) {
        fprintf(stderr, "--tunnel-mole: refusing a multi-conformer atom set\n");
        return 1;
    }
    /* FreeRadius, BRadius and the lining come from the two OPTIONAL trailing
       columns of the atom table (B-factor, atom name). Without the name column
       the extra widths are left equal to the plain radius rather than reported
       as zero, and the header says which happened. */
    mole_pivot_extras(&MA, mbfac, mfree, mbb);
    if (mole_pivot_residues(&MA, mres, &MR) < 0) return 1;
    mole_profile_extras(MA.n ? mbfac : NULL, MA.has_names ? mfree : NULL);
    /* MOLE's DH triangulation by default: same cells as vor_delaunay, plus
       MOLE's per-cell vertex order, which decides lining layer boundaries.
       MOLE_TRIANGULATION=vor selects vor_delaunay, the independent oracle for
       the cell set. */
    if (getenv("MOLE_TRIANGULATION") && !strcmp(getenv("MOLE_TRIANGULATION"), "vor")) {
        if (dt_build(&m, mpiv, np) != 0) { fprintf(stderr, "--tunnel-mole: triangulation failed\n"); return 1; }
        if (mole_build(&M, &m, mpiv, mrad, &P) != 0) { fprintf(stderr, "--tunnel-mole: out of memory\n"); return 1; }
    } else {
        dh_mesh dm;
        if (dh_build(&dm, mpiv, np) != 0) { fprintf(stderr, "--tunnel-mole: triangulation failed\n"); return 1; }
        if (mole_build_from_tetra(&M, dm.tv, dm.tn, dm.nt, mpiv, mrad, &P) != 0) {
            fprintf(stderr, "--tunnel-mole: out of memory\n"); dh_mesh_free(&dm); return 1;
        }
        dh_mesh_free(&dm);
    }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);
    if (getenv("MOLE_TETRA_DUMP")) {
        FILE *tf = fopen("/tmp/ours_tetra.txt", "w");
        int q, kk;
        if (tf) {
            for (q = 0; q < M.nt; q++) {
                int nn = 0;
                for (kk = 0; kk < 4; kk++) if (M.tn[4*q+kk] < 0) nn++;
                fprintf(tf, "%.4f %.4f %.4f %d %d %.6f", M.center[3*q],
                        M.center[3*q+1], M.center[3*q+2], nn, M.boundary[q],
                        M.volume[q]);
                /* The four atoms in vertex order, in MOLE_TETRA_DUMP's field
                   shape, so the dumps compare on order as well as membership. */
                for (kk = 0; kk < 4; kk++) {
                    int a = M.tv[4*q+kk];
                    fprintf(tf, "%c%.3f,%.3f,%.3f", kk ? ';' : ' ',
                            M.axyz[3*a], M.axyz[3*a+1], M.axyz[3*a+2]);
                }
                fputc('\n', tf);
            }
            fclose(tf);
        }
    }
    if (getenv("MOLE_VPROBE")) {
        /* The tetrahedron whose CIRCUMCENTRE is nearest a point, with its four
           atom positions in our vertex order - for the t=0 comparison. */
        double px = 0, py = 0, pz = 0, bd = 1e300; int best = -1, q;
        if (sscanf(getenv("MOLE_VPROBE"), "%lf %lf %lf", &px, &py, &pz) != 3) {
            fprintf(stderr, "MOLE_VPROBE: need three coordinates\n");
        } else {
        for (q = 0; q < M.nt; q++) {
            double dx = M.vcenter[3*q]-px, dy = M.vcenter[3*q+1]-py, dz = M.vcenter[3*q+2]-pz;
            double d2 = dx*dx+dy*dy+dz*dz;
            if (d2 < bd) { bd = d2; best = q; }
        }
        if (best >= 0) {
            fprintf(stderr, "VPROBE idx=%d dist=%.3g vcenter=%.17g %.17g %.17g\n",
                    best, sqrt(bd), M.vcenter[3*best], M.vcenter[3*best+1], M.vcenter[3*best+2]);
            for (q = 0; q < 4; q++) {
                int v = M.tv[4*best+q];
                fprintf(stderr, "  v%d = %.17g %.17g %.17g\n", q,
                        mpiv[3*v], mpiv[3*v+1], mpiv[3*v+2]);
            }
        }
        }
    }
    if (getenv("MOLE_TETRA_PROBE")) {
        /* Everything about the tetrahedron nearest a given centroid, to compare
           against MOLE for the one surface-membership difference. */
        double px = 0, py = 0, pz = 0, bd = 1e300; int best = -1, q;
        if (sscanf(getenv("MOLE_TETRA_PROBE"), "%lf %lf %lf", &px, &py, &pz) != 3)
            fprintf(stderr, "MOLE_TETRA_PROBE: need three coordinates\n");
        else
        for (q = 0; q < M.nt; q++) {
            double dx = M.center[3*q]-px, dy = M.center[3*q+1]-py, dz = M.center[3*q+2]-pz;
            double d2 = dx*dx+dy*dy+dz*dz;
            if (d2 < bd) { bd = d2; best = q; }
        }
        if (best >= 0) {
            double dx = M.center[3*best]-M.vcenter[3*best];
            double dy = M.center[3*best+1]-M.vcenter[3*best+1];
            double dz = M.center[3*best+2]-M.vcenter[3*best+2];
            fprintf(stderr, "PROBE idx=%d dist=%.6f centroid=(%.4f %.4f %.4f)\n",
                    best, sqrt(bd), M.center[3*best], M.center[3*best+1], M.center[3*best+2]);
            fprintf(stderr, "  maxclear=%.17g  2*probe=%.17g  ->%s\n",
                    M.maxclear[best], 2*P.probe_radius,
                    M.maxclear[best] > 2*P.probe_radius ? " clear>2r TRUE" : " clear>2r false");
            fprintf(stderr, "  |c-vc|^2=%.17g  r^2=%.17g  ->%s\n",
                    dx*dx+dy*dy+dz*dz, P.probe_radius*P.probe_radius,
                    (dx*dx+dy*dy+dz*dz) > P.probe_radius*P.probe_radius ? " offset>r TRUE" : " offset>r false");
            fprintf(stderr, "  surface=%d alive=%d boundary=%d comp=%d\n",
                    M.surface[best], M.alive[best], M.boundary[best], M.comp[best]);
            fprintf(stderr, "  neighbours: %d %d %d %d   vertices: %d %d %d %d  (np=%d)\n",
                    M.tn[4*best], M.tn[4*best+1], M.tn[4*best+2], M.tn[4*best+3],
                    M.tv[4*best], M.tv[4*best+1], M.tv[4*best+2], M.tv[4*best+3], np);
            {   int nb2 = 0, kk2;
                for (kk2 = 0; kk2 < 4; kk2++) if (M.tn[4*best+kk2] < 0) nb2++;
                fprintf(stderr, "  missing neighbours: %d\n", nb2);
            }
        }
    }
    if (getenv("MOLE_SURFACE_DEBUG")) {
        /* SurfaceCavity size and its boundary facets, for comparison with
           MOLE's own numbers. A facet is boundary when its neighbour is not a
           surface member; MOLE builds them only for tetrahedra of degree < 4. */
        int nsurf = 0, nfacet = 0, q, kk, iso = 0;
        /* Cavity.Volume = graph.Vertices.Sum(f => f.Volume), so the same
           vertex set the facet count uses. This is the number MOLE prints as
           the MolecularSurface cavity's Volume in cavities.xml, and it is the
           one positive observable that separates a real strict-interior run
           from a crashed one. */
        double svol = 0.0;
        for (q = 0; q < M.nt; q++) {
            int deg = 0, d2;
            if (!M.surface[q]) continue;
            /* AddVerticesAndEdgeRange(cavityGraph.Edges) builds the surface
               graph from EDGES, so a tetrahedron with no surviving neighbour
               never becomes a vertex. */
            for (d2 = 0, kk = 0; kk < 4; kk++) {
                int nb2 = M.tn[4*q+kk];
                if (nb2 >= 0 && M.surface[nb2]) d2++;
            }
            if (d2 == 0) { iso++; continue; }
            nsurf++;
            svol += M.volume[q];
            for (kk = 0; kk < 4; kk++) {
                int nb2 = M.tn[4*q+kk];
                if (nb2 >= 0 && M.surface[nb2]) deg++;
            }
            if (deg >= 4) continue;
            for (kk = 0; kk < 4; kk++) {
                int nb2 = M.tn[4*q+kk];
                if (nb2 < 0 || !M.surface[nb2]) nfacet++;
            }
        }
        fprintf(stderr, "SURFACE tetras=%d boundaryFacets=%d volume=%.3f"
                " (isolated dropped=%d, n_surface=%d)\n",
                nsurf, nfacet, svol, iso, M.n_surface);
        {   /* membership dump, keyed on the centroid, for the MOLE comparison */
            FILE *sf = fopen("/tmp/ours_surface.txt", "w");
            if (sf) {
                for (q = 0; q < M.nt; q++) {
                    int d3 = 0;
                    if (!M.surface[q]) continue;
                    for (kk = 0; kk < 4; kk++) {
                        int nb3 = M.tn[4*q+kk];
                        if (nb3 >= 0 && M.surface[nb3]) d3++;
                    }
                    fprintf(sf, "%.4f %.4f %.4f %d\n", M.center[3*q],
                            M.center[3*q+1], M.center[3*q+2], d3);
                }
                fclose(sf);
            }
        }
    }

    /* ComplexComputation.cs:436-450 numbers CHANNELS by DESCENDING VOLUME and
       assigns Cavity.Id from that; TunnelComparer then orders on the Id. It is
       NOT the order components come out of the pipeline in - on 1ERI the two
       disagree and MOLE's own ids read C1 8.99, C2 7.74, C3 11.94, C4 9.92.
       Stable, so equal volumes keep index order. */
    crank = malloc((size_t)(nc > 0 ? nc : 1) * sizeof(int));
    if (!crank) return 1;
    for (i = 0; i < nc; i++) {
        int r = 0, q;
        for (q = 0; q < nc; q++) {
            if (!cav[q].has_boundary) continue;
            if (cav[q].volume > cav[i].volume
                || (cav[q].volume == cav[i].volume && q < i)) r++;
        }
        crank[i] = r;
    }
    fprintf(stderr, "tunnel_mole: %d atoms, %d tetrahedra, %d channels, %d voids\n",
            MA.n, M.nt, nch, nvd);

    /* GetPaths is its own entry point in MOLE and does not compute tunnels;
       with --path the cavity search is skipped so the output is the paths. */
    for (c = 0; !have_path && c < nc; c++) {
        int origins[MOLE_MAX_ORIGINS], nor, *open = NULL, nop, o, base_;
        double *dist; int *prev;
        if (!(cav[c].has_boundary && cav[c].depth_length > P.min_depth_length
              && cav[c].depth > P.min_depth)) continue;

        if (have_origin) {
            /* A user origin replaces the computed ones for the cavity it lands
               in: nearest vertex by circumcentre with Depth >= 5, rejected past
               OriginRadius - Cavity.GetOrigin's own rule. */
            double bd = 1e300; int best = -1;
            for (i = 0; i < M.nt; i++) {
                double dx, dy, dz, d2;
                if (!M.alive[i] || M.comp[i] != c || M.depth[i] < 5) continue;
                dx = M.vcenter[3*i]-ox; dy = M.vcenter[3*i+1]-oy; dz = M.vcenter[3*i+2]-oz;
                d2 = dx*dx+dy*dy+dz*dz;
                if (d2 < bd) { bd = d2; best = i; }
            }
            if (best < 0 || bd > origin_radius*origin_radius) continue;
            origins[0] = best; nor = 1;
        } else {
            nor = mole_auto_origins(&M, c, auto_cover, max_origins, origins);
            if (getenv("MOLE_ORIGIN_DEBUG"))
                for (i = 0; i < nor; i++)
                    fprintf(stderr, "ORIGIN cav %d depth %d  %.3f %.3f %.3f\n",
                            crank[c] + 1, M.depth[origins[i]],
                            M.center[3*origins[i]], M.center[3*origins[i]+1],
                            M.center[3*origins[i]+2]);   /* Tetrahedron.Center */
        }
        if (!nor) continue;
        /* FilterBoundaryLayers needs Cavity.HasBoundaryVertex, i.e. membership
           in this cavity's BOUNDARY RESIDUE set - the same set cavities.xml
           reports. Built once per cavity rather than per tunnel. */
        if (filter_boundary) {
            int *b = xa_malloc((size_t)(MR.n ? MR.n : 1) * sizeof(int));
            int *n2 = xa_malloc((size_t)(MR.n ? MR.n : 1) * sizeof(int));
            int packed;
            free(bres); bres = xa_calloc((size_t)(MR.n ? MR.n : 1), 1);
            if (b && n2 && bres) {
                int _nb = 0, _ni = 0;
                packed = mole_cavity_residues(&M, c, mres, &MR, b, n2, &_nb, &_ni);
                if (packed >= 0) for (i = 0; i < _nb; i++) bres[b[i]] = 1;
            }
            free(b); free(n2);
        }
        open = NULL; nop = 0;
        if (have_exit && exits_only) {
            /* UseCustomExitsOnly: this cavity's openings are the user exits it
               snapped, and nothing else. */
            int pv = mole_cavity_opening(&M, c, NULL, exitp, origin_radius);
            /* A cavity the exit misses contributes no tunnels of its own, but
               it still yields its origin and the SurfaceCavity is a source from
               it - MOLE prints `SOURCE cavity Id=1 openings=0` followed by the
               surface for the same origin. Skipping the cavity here lost that
               tunnel entirely. */
            if (pv >= 0) {
                open = malloc(sizeof(int));
                if (!open) continue;
                open[0] = pv; nop = 1;
            }
        } else {
            nop = mole_openings(&M, c, surf_cover, &open);
            /* s.Openings holds user exits alongside the computed ones;
               UseCustomExitsOnly filters that list rather than replacing it, so
               a snapped exit is an extra opening. */
            if (have_exit) {
                int pv = mole_cavity_opening(&M, c, NULL, exitp, origin_radius);
                if (pv >= 0) {
                    int *g = xa_realloc(open, (size_t)(nop + 1) * sizeof(int));
                    if (g) { open = g; open[nop++] = pv; }
                }
            }
        }
        if (getenv("MOLE_OPENING_DEBUG")) {
            int z;
            fprintf(stderr, "OPENINGS cav %d n %d :", crank[c] + 1, nop);
            for (z = 0; z < nop; z++) fprintf(stderr, " %d", open[z]);
            fputc('\n', stderr);
        }
        if (!nop && !have_exit) { free(open); continue; }
        dist = malloc((size_t)M.nt*sizeof(double));
        prev = malloc((size_t)M.nt*sizeof(int));
        if (!dist || !prev) return 1;

        for (o = 0; o < nor; o++) {
            base_ = nres;
            mole_dijkstra(&M, c, origins[o], dist, prev, P.weight);
            for (i = 0; i < nop; i++) {
                int v = open[i], len = 0, *p, j, ncp, *cp;
                mole_tunnel_profile pr;
                if (dist[v] >= 1e299) continue;
                for (k = v; k >= 0; k = prev[k]) len++;
                p = malloc((size_t)len*sizeof(int));
                if (!p) continue;
                for (k = v, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
                for (j = 0; j < len; j++) if (M.boundary[p[j]]) break;
                /* Tunnel.Create's FilterTunnelBoundaryLayers: keep a tetrahedron
                   that touches no boundary residue or is deeper than 3 A; of the
                   shallow boundary ones keep at most the FIRST and stop there. */
                if (filter_boundary && bres) {
                    int bc = 0, q, kept = 0;
                    for (q = 0; q < j; q++) {
                        int hasb = 0, z;
                        for (z = 0; z < 4; z++) {
                            int r = mres[M.tv[4*p[q]+z]];
                            if (r >= 0 && bres[r]) { hasb = 1; break; }
                        }
                        if (!hasb || M.depthlen[p[q]] > 3.0) { p[kept++] = p[q]; continue; }
                        if (++bc <= 1) p[kept++] = p[q];
                        else break;
                    }
                    j = kept;
                }
                if (j < 2) { free(p); continue; }
                cp = xa_malloc((size_t)j*sizeof(int));
                ncp = cp ? mole_control_path(&M, p, j, mpiv, mrad, np,
                                             P.interior_threshold, cp) : 0;
                if (ncp < 5) { free(cp); free(p); continue; }
                if (mole_profile(&M, cp, ncp, mpiv, mrad, np, &pr) != 0) {
                    free(cp); free(p); continue;
                }
                if (pr.length >= P.min_tunnel_length
                    && mole_filter_bottleneck(&pr, bottleneck, bottle_tol, 8.0)) {
                    if (nres == ncap) { ncap = ncap ? ncap*2 : 32;
                        res = xa_realloc(res, (size_t)ncap*sizeof(*res)); }
                    res[nres].length = pr.length; res[nres].plen = j;
                    /* Upstream ids are C0 for the SurfaceCavity then C1.. for the
                       ranked cavities, so a regular cavity is crank + 1 and 0 is
                       reserved for the surface. */
                    res[nres].cav = crank[c] + 1;
                    res[nres].npath = j;
                    res[nres].path = xa_malloc((size_t)(j ? j : 1) * sizeof(int));
                    if (res[nres].path) memcpy(res[nres].path, p, (size_t)j * sizeof(int));
                    else res[nres].npath = 0;
                    if (getenv("MOLE_ORIGIN_DEBUG"))
                        fprintf(stderr, "KEEP cav %d origin %d len %.2f\n",
                                crank[c] + 1, o, pr.length);
                    res[nres].prof = pr; nres++;
                } else mole_profile_free(&pr);
                free(cp); free(p);
            }
            /* Per origin, as GetTunnels does - never across origins.
               DEFERRED when a custom exit is in play: GetTunnels builds
               `sources` as this cavity CONCATENATED with the SurfaceCavity and
               runs FilterTunnels over the COMBINED list, so filtering the two
               groups separately keeps a tunnel MOLE removes. The surface branch
               below runs the filter once both are in. */
            if (have_exit)
                surface_tunnel(&M, &P, origins[o], exitp, origin_radius,
                               bottleneck, bottle_tol, &res, &nres, &ncap, np);
            {
                int cnt = nres - base_, q, w;
                if (cnt > 1) {
                    char *dd = xa_calloc((size_t)cnt, 1);
                    mole_tunnel_profile *pf = xa_malloc((size_t)cnt*sizeof(*pf));
                    int *pl = xa_malloc((size_t)cnt*sizeof(int));
                    for (q = 0; q < cnt; q++) { pf[q] = res[base_+q].prof; pl[q] = res[base_+q].plen; }
                    mole_filter_similar(pf, pl, dd, cnt, max_sim);
                    w = base_;
                    for (q = 0; q < cnt; q++) {
                        if (dd[q]) {
                            /* `path` is allocated alongside `prof` and the final
                               cleanup walks SURVIVORS only, so a tunnel dropped
                               here leaked its path array. */
                            mole_profile_free(&res[base_+q].prof);
                            free(res[base_+q].path);
                            res[base_+q].path = NULL; res[base_+q].npath = 0;
                        } else res[w++] = res[base_+q];
                    }
                    nres = w;
                    free(dd); free(pf); free(pl);
                }
            }
        }
        free(dist); free(prev); free(open);
    }
    out = fopen(out_file, "w");
    if (!out) { fprintf(stderr, "--tunnel-mole: cannot write %s\n", out_file); return 1; }
    fprintf(out, "# tunnel-mole %d atoms %d tetrahedra probe %.3f interior %.3f\n",
            MA.n, M.nt, P.probe_radius, P.interior_threshold);
    fprintf(out, "# T id bottleneck length cost throughput ; P id x y z r freeradius bradius\n");
    fprintf(out, "# L id layer start end r freer br localmin charge ionizable npos nneg"
                 " hydropathy hydrophobicity polarity logp logd logs mutability nres"
                 " resn:seq:chain:backbone:flowindex... ; F id nflow resn:seq:chain:backbone...\n");
    fprintf(out, "# H id nhet resn:seq:chain...   (Tunnel.FindHetResidues)\n");
    fprintf(out, "# Y/W id charge ionizable npos nneg hydropathy hydrophobicity polarity"
                 " logp logd logs mutability   (Y unweighted, W weighted by layer length)\n");
    fprintf(out, "# V id type volume depth depthlength nboundary ninner ;"
                 " VB/VI id charge ionizable npos nneg hydropathy hydrophobicity"
                 " polarity logp logd logs mutability nres resn:seq:chain...\n");
    fprintf(out, "# freeradius %s ; bradius %s\n",
            MA.has_names ? "vs backbone+het" : "= radius (no atom-name column in the atom table)",
            "from B-factors (0 where absent)");
    /* Complex.GetPaths: Dijkstra between two user points on the SurfaceCavity
       graph. Both endpoints are given, so no origins and no openings; the
       endpoints snap with GetTetrahedron (nearest member by centroid over the
       WHOLE cavity), not with GetOpening. */
    if (have_path) {
        int ta = mole_cavity_tetrahedron(&M, -1, M.surface, patha, origin_radius);
        int tb = mole_cavity_tetrahedron(&M, -1, M.surface, pathb, origin_radius);
        if (getenv("MOLE_EXIT_DEBUG"))
            fprintf(stderr, "PATH endpoints ta=%d tb=%d\n", ta, tb);
        if (ta >= 0 && tb >= 0) {
            double *dist = xa_malloc((size_t)M.nt*sizeof(double));
            int *prev = xa_malloc((size_t)M.nt*sizeof(int));
            if (dist && prev) {
                mole_dijkstra_mask(&M, -1, ta, dist, prev, P.weight, M.surface);
                if (dist[tb] < 1e299) {
                    int len = 0, *p, j, ncp, *cp, tl, maxd = 0;
                    mole_tunnel_profile pr;
                    for (k = tb; k >= 0; k = prev[k]) len++;
                    p = xa_malloc((size_t)len*sizeof(int));
                    if (p) {
                        for (k = tb, j = len-1; k >= 0; k = prev[k]) p[j--] = k;
                        /* CreatePath trims after the LAST non-boundary
                           tetrahedron - not at the first, as Tunnel.Create
                           does - so interior boundary cells survive mid-route. */
                        for (tl = len-1; tl >= 0; tl--) if (!M.boundary[p[tl]]) break;
                        if (tl >= 10) {
                            len = tl + 1;
                            for (j = 0; j < len; j++)
                                if (M.depth[p[j]] > maxd) maxd = M.depth[p[j]];
                            if (maxd >= 6) {
                                cp = xa_malloc((size_t)len*sizeof(int));
                                ncp = cp ? mole_control_path_ex(&M, p, len, mpiv, mrad, np,
                                                                P.interior_threshold, cp, 1) : 0;
                                if (getenv("MOLE_EXIT_DEBUG"))
                                    fprintf(stderr, "PATH path=%d controlPath=%d maxDepth=%d\n",
                                            len, ncp, maxd);
                                if (ncp >= 5 && mole_profile(&M, cp, ncp, mpiv, mrad, np, &pr) == 0) {
                                    if (nres == ncap) { ncap = ncap ? ncap*2 : 32;
                                        res = xa_realloc(res, (size_t)ncap*sizeof(*res)); }
                                    res[nres].length = pr.length; res[nres].plen = len;
                                    res[nres].cav = 0; res[nres].prof = pr;
                                    res[nres].npath = len;
                                    res[nres].path = xa_malloc((size_t)(len ? len : 1) * sizeof(int));
                                    if (res[nres].path) memcpy(res[nres].path, p, (size_t)len * sizeof(int));
                                    else res[nres].npath = 0;
                                    nres++;
                                }
                                free(cp);
                            }
                        }
                        free(p);
                    }
                }
            }
            free(dist); free(prev);
        }
    }

    /* The SurfaceCavity. GetTunnels builds `sources` as the origin's cavity
       concatenated with it on every run - UseCustomExitsOnly only filters each
       source's openings. Listed after the regular cavities
       (TunnelComputation.cs:91)
       and it has no openings of its own - UpdateOpenings returns immediately
       for it - so it yields tunnels only through a user exit. */
    /* the SurfaceCavity is a source per ORIGIN; see surface_tunnel */



    /* TunnelCollection.TunnelComparer, applied once every branch has appended.
       The SurfaceCavity is Id 0, so its tunnels sort first. */
    if (nres > 1) qsort(res, (size_t)nres, sizeof(*res), by_cav_len);

    /* Cavities first: they describe the space the tunnels were found in, and
       MOLE exports them as their own file. Only the ones that passed the depth
       filters, which are the ones tunnels could come from. */
    if (MA.has_names) {
        for (c = 0; c < nc; c++) {
            if (!(cav[c].has_boundary && cav[c].depth_length > P.min_depth_length
                  && cav[c].depth > P.min_depth)) continue;
            /* crank is the descending-volume rank MOLE assigns Cavity.Id from,
               so the ids line up with its cavities.xml rather than being a
               sequential counter over whatever survived the depth filters. */
            write_cavity(out, crank[c] + 1, &M, c, &cav[c], mres);
        }
    }
    for (i = 0; i < nres; i++) {
        int ns = (int)(res[i].length * 8.0), q;
        double dt, bmin = 1e300;
        if (ns < 1) ns = 1;
        dt = 1.0/ns;
        for (q = 0; q <= ns; q++) {
            double r = mole_spline_eval(&res[i].prof.sr, dt*q);
            if (r < bmin) bmin = r;
        }
        /* Field 6 is the CAVITY id, upstream's own numbering: 0 is the
           SurfaceCavity (upstream "C0") and 1.. are the ranked buried cavities,
           so a tunnel can be labelled T<n>C<cav> the way MOLE itself does. It was
           tracked internally all along and simply never written, leaving
           auto-origin runs unable to say WHICH cavity produced a tunnel or to
           distinguish a surface-cavity route. Appended, not substituted: every
           existing reader indexes fields 2-5 positionally. Columns 4 and 5 stay
           0.0 - MOLE has no cost/throughput concept (see _tunnel_throughput,
           which derives throughput from the profile instead). */
        fprintf(out, "T %d %.4f %.4f %.6f %.6f %d\n",
                i+1, bmin, res[i].length, 0.0, 0.0, res[i].cav);
        for (q = 0; q <= ns; q++) {
            double t = dt*q;
            fprintf(out, "P %d %.4f %.4f %.4f %.4f %.4f %.4f\n", i+1,
                    mole_spline_eval(&res[i].prof.sx, t),
                    mole_spline_eval(&res[i].prof.sy, t),
                    mole_spline_eval(&res[i].prof.sz, t),
                    mole_spline_eval(&res[i].prof.sr, t),
                    mole_spline_eval(&res[i].prof.sfr, t),
                    mole_spline_eval(&res[i].prof.sbr, t));
        }
        if (MA.has_names) write_lining(out, i + 1, &res[i].prof, np);
        write_het(out, i + 1, &res[i].prof, &M, res[i].path, res[i].npath, np);
    }
    fclose(out);
    fprintf(stderr, "tunnel_mole: %d tunnels\n", nres);
    for (i = 0; i < nres; i++) mole_profile_free(&res[i].prof);
    for (i = 0; i < nres; i++) free(res[i].path);
    free(res); free(crank); free(bres); free(cav); mole_free(&M); dt_free(&m);
    return 0;
}
