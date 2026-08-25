/* MOLE 2's cavity pipeline: atoms -> Voronoi graph -> cavity diagrams.
 *
 * Ported from WebChemistry.Tunnels.Core (MIT), following the SOURCE rather than
 * the two method papers, which disagree with it in three places recorded in
 * Every quantity here has a named counterpart in MOLE 2's own source;
 * where the arithmetic form matters for reproducing their numbers it is written
 * the way they write it rather than the way it would otherwise be natural.
 *
 * Checkpoint for this stage: 9 channels and 6 voids on 1tqn.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "mole_complex.h"
#include "mole_rng.h"
#include "../xalloc.h"

/* MOLE's own table (TunnelVdwRadii.cs). NOT AMBER, despite both papers saying
   so. Anything absent falls back to a general element table; 1tqn uses only
   C/N/O/S/FE, all of which are listed here. */
double mole_vdw_radius(const char *elem)
{
    /* MOLE's tunnel-specific overrides (TunnelVdwRadii.cs). NOT AMBER, despite
       both papers saying so. */
    static const struct { const char *e; double r; } TUN[] = {
        { "H", 1.0 }, { "O", 1.45 }, { "S", 1.77 }, { "N", 1.55 }, { "C", 1.61 },
        { "FE", 1.7 }, { "P", 1.7 }, { "SI", 1.8 }, { "AL", 1.84 },
        { "LI", 1.8 }, { "NA", 1.8 }, { "CL", 1.75 }
    };
    /* Anything else falls through to the general element table
       (ElementAndBondInfo, populated in BondInfo.cs), NOT to a flat constant.
       Returning one value for every unlisted element was wrong for any
       selection containing Zn, Ca, Mg, Br and so on - which is the normal case
       once a user picks their own atoms in VMD. Note P, Li and Na differ
       between the two tables (1.7/1.9, 1.8/1.82, 1.8/2.27); the tunnel table
       wins, which is why it is consulted first. */
    static const struct { const char *e; double r; } GEN[] = {
    { "AS", 1.85 },
    { "B", 0.9 },
    { "BE", 1.53 },
    { "BR", 1.85 },
    { "CA", 2.31 },
    { "D", 1.0 },
    { "HF", 1.59 },
    { "HG", 1.0 },
    { "I", 1.85 },
    { "IR", 1.47 },
    { "K", 2.75 },
    { "MG", 1.73 },
    { "MN", 1.0 },
    { "PT", 1.0 },
    { "RH", 1.49 },
    { "RU", 1.34 },
    { "SE", 1.9 },
    { "SR", 2.49 },
    { "TE", 1.0 },
    { "W", 1.39 }
    };
    size_t i;
    char u[8];
    for (i = 0; i + 1 < sizeof u && elem[i]; i++)
        u[i] = (elem[i] >= 'a' && elem[i] <= 'z') ? (char)(elem[i] - 32) : elem[i];
    u[i] = 0;
    for (i = 0; i < sizeof TUN / sizeof TUN[0]; i++)
        if (!strcmp(u, TUN[i].e)) return TUN[i].r;
    for (i = 0; i < sizeof GEN / sizeof GEN[0]; i++)
        if (!strcmp(u, GEN[i].e)) return GEN[i].r;
    return 1.75;   /* ElementInfo.Default.VdwRadius */
}

/* MOLE's own decimal parser (NumberParser.ParseDoubleFast), reproduced.
 *
 * It accumulates the integer part, then the fractional digits as an integer over
 * a power of ten, and returns main + point/div. That rounds TWICE - once on the
 * division, once on the addition - where a correctly-rounded strtod rounds once.
 * On 1tqn it puts 252 of 11997 coordinates one ULP away from the correctly
 * rounded value, e.g. "6.919" becomes 6.9190000000000005.
 *
 * This is upstream of everything: their triangulation is built on these numbers,
 * so a port that parses correctly is already working from different input. Only
 * MOLE mode needs it - in the plugin proper VMD supplies the coordinates. */
double mole_parse_double(const char *s)
{
    double main_ = 0.0, point = 0.0, div = 1.0;
    int neg = 0;
    const char *p = s;
    for (; *p; p++) {
        if (*p >= '0' && *p <= '9') main_ = main_ * 10 + (*p - '0');
        else if (*p == '-') neg = 1;
        else if (*p == '.') break;
        else if (*p == 'e' || *p == 'E') return neg ? -main_ : main_;
        else if (*p == ' ' || *p == '\t') { if (main_ != 0.0) break; }
        else break;
    }
    if (*p == '.') {
        for (p++; *p >= '0' && *p <= '9'; p++) { div *= 10; point = point * 10 + (*p - '0'); }
    }
    return neg ? -(main_ + point / div) : (main_ + point / div);
}

/* Case-insensitive membership in a NULL-terminated name table. */
static int name_in(const char *s, const char *const *tab)
{
    size_t i;
    char u[8];
    for (i = 0; i + 1 < sizeof u && s[i]; i++)
        u[i] = (s[i] >= 'a' && s[i] <= 'z') ? (char)(s[i] - 32) : s[i];
    u[i] = 0;
    for (i = 0; tab[i]; i++) if (!strcmp(u, tab[i])) return 1;
    return 0;
}

/* PdbEx.backboneNames - amino AND nucleic, and note it includes H. This is a
   NAME test, not an element test: "CA" here is the alpha carbon, while the very
   same string in the vdW table is calcium. */
int mole_is_backbone_name(const char *name)
{
    static const char *const BB[] = {
        "C", "N", "O", "H", "CA", "P", "O1P", "O2P", "OP1", "OP2",
        "O5'", "C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O3'", "O2'", NULL };
    return name_in(name, BB);
}

/* PdbResidue.aminoNames. */
int mole_is_amino_name(const char *resn)
{
    static const char *const AA[] = {
        "ALA", "ARG", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
        "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ASN",
        NULL };
    return name_in(resn, AA);
}

/* Read the atom table the oracle comparison uses:
   "x y z element is_water chain seq resname", one atom per line, FILE ORDER.
   Coordinates go through mole_parse_double so the input is bit-identical to
   MOLE's, which its own parser makes non-obvious. In the plugin proper VMD
   supplies these and none of this is needed. */
/* Read one frame's atoms from VMDHole's packed record instead of the 12-column
   text. Same atoms, same doubles: the writer sends raw IEEE coordinates and the
   quantisation happens HERE, by rendering %.3f and running it back through the
   very parser the text path uses. That is what keeps the two inputs bit-identical
   - MOLE's parser is a hand-rolled main+point/div, not strtod, so handing it raw
   doubles would shift coordinates by an ULP and Delaunay predicates are decided
   at exactly that scale. Verified on 56031 real coordinates that this C %.3f and
   the writer's Tcl %.3f agree character for character.
   Layout, little-endian, written and read on the same host:
     "MOLEATM1" | int32 version | int32 natoms | int32 has_names
     natoms x (elem[8] chain[8] resn[8] name[8] alt[4] icode[4])
     natoms x int32 water | natoms x int32 seq | natoms x double bfac
     natoms x double x | natoms x double y | natoms x double z          */
static int mole_read_atoms_packed(FILE *f, mole_atoms *A, int cap,
                                  const char *path, int exact)
{
    int ver = 0, n = 0, hn = 0, i;
    char buf[64];
    if (fread(&ver, 4, 1, f) != 1 || fread(&n, 4, 1, f) != 1
        || fread(&hn, 4, 1, f) != 1) return -1;
    if (ver != 1) {
        fprintf(stderr, "mole_read_atoms: %s: unsupported packed version %d\n", path, ver);
        return -1;
    }
    if (n < 0) return -1;
    if (n > cap) {
        fprintf(stderr, "mole_read_atoms: %s holds more than %d atoms - "
                "refusing to analyse a truncated structure. Narrow the "
                "selection or rebuild with a larger MOLE_MAXATOM.\n", path, cap);
        return -1;
    }
    for (i = 0; i < n; i++) {
        if (fread(A->elem[i], 8, 1, f) != 1) return -1;
        if (fread(A->chain[i], 8, 1, f) != 1) return -1;
        if (fread(A->resn[i], 8, 1, f) != 1) return -1;
        if (fread(A->name[i], 8, 1, f) != 1) return -1;
        if (fread(A->alt[i], 4, 1, f) != 1) return -1;
        if (fread(A->icode[i], 4, 1, f) != 1) return -1;
        A->elem[i][7] = A->chain[i][7] = A->resn[i][7] = A->name[i][7] = 0;
        A->alt[i][3] = A->icode[i][3] = 0;
        /* Same refusal the text reader makes: a bare 0/1 in the name slot is the
           old free flag, and accepting it would leave every atom non-backbone -
           changing FreeRadius and the per-layer backbone split with no error. */
        if (A->name[i][1] == 0 && (A->name[i][0] == '0' || A->name[i][0] == '1')) {
            fprintf(stderr, "mole_read_atoms: %s atom %d: name is a bare 0/1 - "
                    "that was the old free flag; it is now the atom name\n", path, i + 1);
            return -1;
        }
    }
    for (i = 0; i < n; i++) if (fread(&A->water[i], 4, 1, f) != 1) return -1;
    for (i = 0; i < n; i++) if (fread(&A->seq[i],   4, 1, f) != 1) return -1;
    for (i = 0; i < n; i++) if (fread(&A->bfac[i],  8, 1, f) != 1) return -1;
    for (i = 0; i < 3; i++) {
        int j;
        for (j = 0; j < n; j++) {
            double v;
            if (fread(&v, 8, 1, f) != 1) return -1;
            snprintf(buf, sizeof buf, "%.3f", v);
            A->xyz[3*j+i] = exact ? strtod(buf, NULL) : mole_parse_double(buf);
        }
    }
    A->has_names = hn ? 1 : 0;
    for (i = 0; i < n; i++) {
        if (A->has_names) {
            A->backbone[i] = mole_is_backbone_name(A->name[i]);
            A->freeatom[i] = A->backbone[i] || !mole_is_amino_name(A->resn[i]);
        } else {
            A->backbone[i] = 0;
            A->freeatom[i] = 0;
        }
    }
    A->n = n;
    return n;
}

int mole_read_atoms(const char *path, mole_atoms *A, int cap)
{
    FILE *f = fopen(path, "r");
    char sx[32], sy[32], sz[32], line[512];
    int n = 0, exact = getenv("MOLE_ATOMS_EXACT") != NULL;
    char magic[8];
    if (!f) return -1;
    A->has_names = 0;
    /* Keyed on the file's own magic rather than its name, because every caller
       (the engine and six dump/test programs) reaches this one reader. */
    if (fread(magic, 1, 8, f) == 8 && memcmp(magic, "MOLEATM1", 8) == 0) {
        n = mole_read_atoms_packed(f, A, cap, path, exact);
        fclose(f);
        return n;
    }
    rewind(f);
    n = 0;
    /* Line-based rather than a bare fscanf: the two OPTIONAL trailing columns
       (bfactor, atom name) must not desynchronise the parse when absent, which
       a field-counting fscanf across newlines would do. */
    while (n < cap && fgets(line, sizeof line, f)) {
        int nf;
        char sb[32];
        A->bfac[n] = 0.0; A->name[n][0] = 0;
        A->backbone[n] = 0; A->freeatom[n] = 0;
        strcpy(A->alt[n], "-"); strcpy(A->icode[n], "-");
        nf = sscanf(line, "%31s %31s %31s %7s %d %7s %d %7s %31s %7s %3s %3s",
                    sx, sy, sz, A->elem[n], &A->water[n],
                    A->chain[n], &A->seq[n], A->resn[n], sb, A->name[n],
                    A->alt[n], A->icode[n]);
        if (nf < 8) continue;
        if (nf >= 9)  A->bfac[n] = atof(sb);
        if (nf >= 10) {
            /* Column 10 is the atom name. A bare 0/1 is a superseded free flag: a
               stale table would classify every atom non-backbone and change
               FreeRadius, so refuse it. */
            if (A->name[n][1] == 0 && (A->name[n][0] == '0' || A->name[n][0] == '1')) {
                fprintf(stderr, "mole_read_atoms: %s line %d: column 10 is a bare "
                        "0/1 - that was the old free flag; it is now the atom name\n",
                        path, n + 1);
                fclose(f); return -1;
            }
            A->backbone[n] = mole_is_backbone_name(A->name[n]);
            /* PdbEx.IsHetAtom is "HETATM record OR not a standard residue name".
               An atom table has no record type, so only the second half is
               testable. */
            A->freeatom[n] = A->backbone[n] || !mole_is_amino_name(A->resn[n]);
            A->has_names = 1;
        }
        /* MOLE_ATOMS_EXACT: the table came from MOLE's own G17 dump, so it
           already holds their parsed doubles and must be read back with a
           correctly-rounded strtod. Applying their parser again would round a
           second time. Normal .cif-derived tables go through their parser. */
        if (exact) {
            A->xyz[3*n+0] = strtod(sx, NULL);
            A->xyz[3*n+1] = strtod(sy, NULL);
            A->xyz[3*n+2] = strtod(sz, NULL);
        } else {
            A->xyz[3*n+0] = mole_parse_double(sx);
            A->xyz[3*n+1] = mole_parse_double(sy);
            A->xyz[3*n+2] = mole_parse_double(sz);
        }
        n++;
    }
    /* REFUSE a structure that does not fit, rather than analysing its first
       `cap` atoms. Truncating is worse than failing: the run still produces a
       plausible tunnel set, from a partial assembly, with nothing to tell the
       user that the rest of the structure was discarded - so bottlenecks and
       lining residues can both belong to a molecule that was never complete.
       Large biological assemblies are ordinary MOLE input. */
    if (n == cap) {
        char probe[512];
        while (fgets(probe, sizeof probe, f)) {
            char a[32], b[32], c[32], d[8];
            if (sscanf(probe, "%31s %31s %31s %7s", a, b, c, d) >= 4) {
                fprintf(stderr, "mole_read_atoms: %s holds more than %d atoms - "
                        "refusing to analyse a truncated structure. Narrow the "
                        "selection or rebuild with a larger MOLE_MAXATOM.\n",
                        path, cap);
                fclose(f); return -1;
            }
        }
    }
    fclose(f);
    A->n = n;
    return n;
}

/* Per-pivot bfactor, free flag and backbone flag, in the same order
   mole_pivots produced. */
int mole_pivot_extras(const mole_atoms *A, double *bfac, int *freeatom,
                      int *backbone)
{
    int p;
    /* PIVOT order, not file order: mole_pivots reorders, and these arrays are
       indexed by pivot. Requires mole_pivots to have run. */
    for (p = 0; p < A->npiv; p++) {
        int i = A->order[p];
        if (bfac) bfac[p] = A->bfac[i];
        /* No name column at all: every atom counts, so FreeRadius == Radius and
           the caller can say so rather than reporting a radius measured against
           an empty set. */
        if (freeatom) freeatom[p] = A->has_names ? A->freeatom[i] : 1;
        if (backbone) backbone[p] = A->backbone[i];
    }
    return A->npiv;
}

/* Jitter every atom, then keep the non-water ones as triangulation pivots.
   pivsrc maps a pivot back to its atom, which the origin report needs. */
/* MOLE's pivot order: chains ordered by Identifier, then residues by Number/
 * InsertionCode, atoms by alt-loc (blank first) then file order - LINQ's
 * OrderBy, which is stable. Feeds the incremental triangulation, so it fixes
 * the per-cell vertex order mole_dh.c reproduces.
 *
 * Chain identifiers COLLATE under InvariantCulture, not strcmp: case-
 * insensitive first, lowercase before uppercase as a tiebreak ("a" sorts
 * before "B", after "A"). Only the alphanumeric range is reproduced; two
 * DIFFERENT punctuation chain ids (reachable - extract_atoms_pdb.py writes
 * '?' for a blank chain) collate by an ICU weight table and are a declared
 * limit, not an unreachable case.
 */
static int chain_class(unsigned char c)
{
    if (isalpha(c)) return 2;
    if (isdigit(c)) return 1;
    return 0;
}

int mole_chain_cmp(const char *a, const char *b)
{
    int tie = 0;
    for (; *a && *b; a++, b++) {
        unsigned char ua = (unsigned char)*a, ub = (unsigned char)*b;
        int ka = chain_class(ua), kb = chain_class(ub), ca, cb;
        if (ka != kb) return ka < kb ? -1 : 1;
        ca = tolower(ua); cb = tolower(ub);
        if (ca != cb) return ca < cb ? -1 : 1;
        if (!tie && ua != ub) tie = islower(ua) ? -1 : 1;
    }
    if (*a) return 1;
    if (*b) return -1;
    return tie;
}

/* Residue identity as MOLE keys it: chain, number, insertion code. The residue
   NAME is deliberately absent - MOLE merges two names sharing an identifier into
   one residue, and mole_pivot_residues (which does include the name, for the
   lining) is a separate question from ordering. */
typedef struct { int id, seq, first; const char *icode; } piv_res;

static const mole_atoms *g_sort_atoms;
static int piv_res_cmp(const void *pa, const void *pb)
{
    const piv_res *a = pa, *b = pb;
    const mole_atoms *A = g_sort_atoms;
    int c = mole_chain_cmp(A->chain[a->first], A->chain[b->first]);
    if (c) return c;
    if (a->seq != b->seq) return a->seq < b->seq ? -1 : 1;
    c = strcmp(a->icode, b->icode);
    if (c) return c;
    return a->id < b->id ? -1 : 1;     /* LINQ's OrderBy is stable */
}

/* Fills A->order / A->npiv. Returns the pivot count. */
static int mole_pivot_order(mole_atoms *A)
{
    static piv_res R[MOLE_MAXATOM];
    static int rof[MOLE_MAXATOM];      /* residue index per atom, -1 for water */
    static int roff[MOLE_MAXATOM + 1]; /* residue -> slice of ratoms */
    static int ratoms[MOLE_MAXATOM];   /* atom indices, grouped by residue, file order */
    int i, j, nr = 0, np = 0;

    for (i = 0; i < A->n; i++) {
        rof[i] = -1;
        if (A->water[i]) continue;
        /* A residue's atoms are contiguous in every real file, so probe the last
           one first; the scan behind it keeps this correct when they are not. */
        j = nr - 1;
        if (j < 0 || R[j].seq != A->seq[i] ||
            strcmp(A->chain[R[j].first], A->chain[i]) ||
            strcmp(R[j].icode, A->icode[i])) {
            for (j = 0; j < nr; j++)
                if (R[j].seq == A->seq[i] &&
                    !strcmp(A->chain[R[j].first], A->chain[i]) &&
                    !strcmp(R[j].icode, A->icode[i])) break;
        }
        if (j == nr) {
            R[nr].seq = A->seq[i]; R[nr].first = i;
            R[nr].icode = A->icode[i]; R[nr].id = nr;
            nr++;
        } else if (strcmp(A->resn[R[j].first], A->resn[i])) {
            /* MOLE's residue identity has no NAME in it, so it merges these into
               one residue while mole_pivot_residues - which does key on the name,
               for the lining - splits them. Nothing here can resolve that
               disagreement; say so rather than produce a quietly different
               lining. The altLoc warning below is the same guarantee. */
            fprintf(stderr, "mole_pivots: %s %d holds both %s and %s; MOLE treats "
                    "one residue where the lining will see two\n",
                    A->chain[i], A->seq[i], A->resn[R[j].first], A->resn[i]);
        }
        rof[i] = j;
    }

    /* Bucket the atoms by residue in one pass. A scan per residue would be
       O(residues * atoms), which is minutes on a large structure. */
    for (j = 0; j <= nr; j++) roff[j] = 0;
    for (i = 0; i < A->n; i++) if (rof[i] >= 0) roff[rof[i] + 1]++;
    for (j = 0; j < nr; j++) roff[j + 1] += roff[j];
    {
        static int fill[MOLE_MAXATOM];
        for (j = 0; j < nr; j++) fill[j] = roff[j];
        for (i = 0; i < A->n; i++) if (rof[i] >= 0) ratoms[fill[rof[i]]++] = i;
    }

    g_sort_atoms = A;
    qsort(R, (size_t)nr, sizeof R[0], piv_res_cmp);

    for (j = 0; j < nr; j++) {
        int r = R[j].id, lo = roff[r], hi = roff[r + 1], nalt = 0;
        /* Atoms of this residue grouped by altLoc: blank first, then ascending,
           file order inside a group. Repeated passes rather than a sort, so the
           within-group order is the file's without relying on qsort being
           stable, which it is not. */
        const char *want = "-";
        for (;;) {
            const char *next = NULL;
            for (i = lo; i < hi; i++) {
                int a = ratoms[i];
                if (!strcmp(A->alt[a], want)) A->order[np++] = a;
                else if (strcmp(A->alt[a], want) > 0 &&
                         (!next || strcmp(A->alt[a], next) < 0)) next = A->alt[a];
            }
            if (!next) break;
            want = next;
            if (++nalt == 2) {
                /* Two mutually exclusive conformers in one atom set is not a
                   fidelity question, it is a structure that cannot exist: the
                   triangulation would be of a hybrid, changing cavities,
                   bottlenecks and lining. MOLE keeps the blank group plus the
                   first non-blank one; rather than silently pick for the user,
                   refuse and say which selection to make. */
                fprintf(stderr, "mole_pivots: residue %s %d carries more than "
                        "one alternate location. Select one conformer - e.g. "
                        "`altloc \"\" or altloc A` - and re-run.\n",
                        A->chain[R[j].first], R[j].seq);
                A->npiv = -1;
                return -1;
            }
        }
    }
    A->npiv = np;
    return np;
}

int mole_is_het_resname(const char *resn) { return !mole_is_amino_name(resn); }

/* Returns the pivot count, or -1 when the atom set holds more than one
   alternate conformer for a residue - see mole_pivot_order. */
int mole_pivots(mole_atoms *A, double *piv, double *rad, int *pivsrc)
{
    int p, np;
    /* Snapshot BEFORE the jitter: this is MOLE's InvariantPosition, which
       FindHetResidues - and only FindHetResidues - searches. */
    static double raw[3*MOLE_MAXATOM];
    memcpy(raw, A->xyz, (size_t)A->n * 3 * sizeof(double));
    mole_jitter(A->xyz, A->n, 0);      /* drawn in FILE order, before any of this */
    np = mole_pivot_order(A);
    if (np < 0) return -1;
    for (p = 0; p < np; p++) {
        int i = A->order[p];
        memcpy(&A->invpiv[3*p], &raw[3*i], 3*sizeof(double));
        memcpy(&piv[3*p], &A->xyz[3*i], 3*sizeof(double));
        rad[p] = mole_vdw_radius(A->elem[i]);
        if (pivsrc) pivsrc[p] = i;
    }
    return np;
}

double mole_det4(const double m[4][4])
{
    /* Matrix3D.Determinant's exact evaluation - first-row expansion with six
       shared 2x2 minors from the bottom two rows, in this order.
       Matrix3D's field layout maps as
           m11 m12 m13 m14        m[0][0..3]
           m21 m22 m23 m24        m[1][0..3]
           m31 m32 m33 m34        m[2][0..3]
           offX offY offZ m44     m[3][0..3]
       This is a FIDELITY match, not a behaviour fix: measured over 20000 random
       4-atom sets from 1BL8, this and a last-column expansion give bit-identical
       circumcentres every time, so nothing observable depended on it. It is
       written MOLE's way so that a future divergence cannot come from here.
       1BL8's four spurious layer boundaries are caused by per-cell VERTEX
       ORDER, not by this. */
    double num1 = m[2][2]*m[3][3] - m[2][3]*m[3][2];
    double num2 = m[2][1]*m[3][3] - m[2][3]*m[3][1];
    double num3 = m[2][0]*m[3][3] - m[2][3]*m[3][0];
    double num4 = m[2][1]*m[3][2] - m[2][2]*m[3][1];
    double num5 = m[2][0]*m[3][2] - m[2][2]*m[3][0];
    double num6 = m[2][0]*m[3][1] - m[2][1]*m[3][0];

    return m[0][0] * (m[1][1]*num1 - m[1][2]*num2 + m[1][3]*num4) -
           m[0][1] * (m[1][0]*num1 - m[1][2]*num3 + m[1][3]*num5) +
           m[0][2] * (m[1][0]*num2 - m[1][1]*num3 + m[1][3]*num6) -
           m[0][3] * (m[1][0]*num4 - m[1][1]*num5 + m[1][2]*num6);
}

/* MathHelper.SphereFromPoints: the circumcentre as a ratio of 4x4
   determinants. Only the centre is used downstream. */
void mole_circumsphere(const double p[4][3], double *cx, double *cy, double *cz)
{
    double m[4][4], a, dx, dy, dz, s;
    int i;
    double ls[4];
    for (i = 0; i < 4; i++) ls[i] = p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2];

    for (i = 0; i < 4; i++) { m[i][0]=p[i][0]; m[i][1]=p[i][1]; m[i][2]=p[i][2]; m[i][3]=1; }
    a = mole_det4(m);
    for (i = 0; i < 4; i++) m[i][0] = ls[i];
    dx = mole_det4(m);
    for (i = 0; i < 4; i++) m[i][1] = p[i][0];
    dy = -mole_det4(m);
    for (i = 0; i < 4; i++) m[i][2] = p[i][1];
    dz = mole_det4(m);

    s = 1.0 / (2.0 * a);
    *cx = s * dx; *cy = s * dy; *cz = s * dz;
}

static double dist3(const double *a, const double *b)
{
    double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* Vector3D.Normalize: multiply by 1/sqrt(L2). NOT a division by sqrt(L2) -
   the two round differently, and 290 of 766 edge clearances came out 1 ULP
   away when this was written the natural way. */
static void normalize3(double *v)
{
    double inv = 1.0 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] *= inv; v[1] *= inv; v[2] *= inv;
}

/* Plane3D.FromPoints then DistanceTo, which is ABSOLUTE in their code.
   Normal = normalize((b-a) x (c-b)), D = -(c . n). */
static double plane_dist(const double *a, const double *b, const double *c,
                         const double *pt)
{
    double u[3], v[3], n[3], d;
    int i;
    for (i = 0; i < 3; i++) { u[i] = b[i]-a[i]; v[i] = c[i]-b[i]; }
    n[0] = u[1]*v[2] - u[2]*v[1];
    n[1] = u[2]*v[0] - u[0]*v[2];
    n[2] = u[0]*v[1] - u[1]*v[0];
    normalize3(n);
    d = c[0]*n[0] + c[1]*n[1] + c[2]*n[2];
    return fabs(n[0]*pt[0] + n[1]*pt[1] + n[2]*pt[2] - d);
}

/* Distance from a point to the INFINITE line through a and b
   (MathEx.DistanceTo(this Line3D, ...) - the direction is normalised, so this
   is not a segment distance). */
static double line_dist(const double *a, const double *b, const double *pt)
{
    double d[3], w[3], c[3];
    int i;
    for (i = 0; i < 3; i++) { d[i] = b[i]-a[i]; w[i] = a[i]-pt[i]; }
    normalize3(d);
    c[0] = d[1]*w[2] - d[2]*w[1];
    c[1] = d[2]*w[0] - d[0]*w[2];
    c[2] = d[0]*w[1] - d[1]*w[0];
    /* Divided by |d| even though d was just normalised: |normalize(v)| is not
       exactly 1, and they do the division, so it must be done here too. */
    return sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
         / sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
}

/* Tetrahedron.PartialVolume: the vdW sphere's share of a corner, approximated
   by the cone on the three normalised edge directions. */
static double partial_volume(const double p[4][3], const double *r, int i)
{
    double a[3], b[3], c[3], la = 0, lb = 0, lc = 0, cr[3], dp;
    int k;
    for (k = 0; k < 3; k++) {
        a[k] = p[(i+1)%4][k] - p[i][k];
        b[k] = p[(i+2)%4][k] - p[i][k];
        c[k] = p[(i+3)%4][k] - p[i][k];
    }
    normalize3(a); normalize3(b); normalize3(c);
    (void)la; (void)lb; (void)lc;
    cr[0] = b[1]*c[2] - b[2]*c[1];
    cr[1] = b[2]*c[0] - b[0]*c[2];
    cr[2] = b[0]*c[1] - b[1]*c[0];
    dp = a[0]*cr[0] + a[1]*cr[1] + a[2]*cr[2];
    return r[i]*r[i]*r[i] * fabs(dp) / 6.0;
}

/* Tetrahedron.Init + Update(strict, cavityGraph).
 *
 * strict:false - MaxClearance is a MAXIMUM over the six pairwise atom gaps and
 * the four plane-to-opposite-vertex gaps, not a minimum and not a distance to
 * the nearest atom. A tetrahedron counts as interior only when EVERY one of
 * those ten is small.
 *
 * strict:true - MaxClearance is instead the maximum over the ADJACENT EDGES'
 * Clearance, and those clearances are all still zero at this point. MOLE's own
 * comment at ComplexComputation.cs:344 says "if strict, need to update the
 * edges 1st!" and the line that would do it (345/348) is commented out, so
 * Edge.Clearance is the C# double default. The first live Edge.Update is at
 * :356, AFTER RemoveOutsideAndInterior. Reproduced rather than repaired: fixing
 * it here would make the port disagree with the reference tool. */
static void tet_init(mole_complex *M, int t, int strict)
{
    double p[4][3], r[4], d = -1000.0, td;
    int i, j, k;
    for (i = 0; i < 4; i++) {
        int a = M->tv[4*t+i];
        for (k = 0; k < 3; k++) p[i][k] = M->axyz[3*a+k];
        r[i] = M->arad[a];
    }
    for (k = 0; k < 3; k++)
        M->center[3*t+k] = 0.25 * (p[0][k] + p[1][k] + p[2][k] + p[3][k]);
    mole_circumsphere(p, &M->vcenter[3*t], &M->vcenter[3*t+1], &M->vcenter[3*t+2]);

    {
        double total, cr[3], u[3], v[3], w[3];
        for (k = 0; k < 3; k++) { u[k]=p[1][k]-p[0][k]; v[k]=p[2][k]-p[0][k]; w[k]=p[3][k]-p[0][k]; }
        cr[0] = v[1]*w[2] - v[2]*w[1];
        cr[1] = v[2]*w[0] - v[0]*w[2];
        cr[2] = v[0]*w[1] - v[1]*w[0];
        total = fabs(u[0]*cr[0] + u[1]*cr[1] + u[2]*cr[2]) / 6.0;
        for (i = 0; i < 4; i++) total -= partial_volume(p, r, i);
        M->volume[t] = total;
    }

    if (strict) {
        /* graph.AdjacentEdgeList(this): one edge per FINITE neighbour. A
           tetrahedron with none is not a graph vertex at all and MOLE throws
           there; d keeps its -1000 sentinel instead, which the interior test
           reads the same way as the zero it would otherwise see. */
        for (k = 0; k < 4; k++) {
            if (M->tn[4*t+k] < 0) continue;
            td = M->eclear[4*t+k];
            if (td > d) d = td;
        }
    } else {
        for (i = 0; i < 3; i++)
            for (j = i+1; j < 4; j++) {
                td = dist3(p[i], p[j]) - r[i] - r[j];
                if (td > d) d = td;
            }
        td = plane_dist(p[0], p[1], p[2], p[3]) - r[3]; if (td > d) d = td;
        td = plane_dist(p[0], p[1], p[3], p[2]) - r[2]; if (td > d) d = td;
        td = plane_dist(p[0], p[2], p[3], p[1]) - r[1]; if (td > d) d = td;
        td = plane_dist(p[1], p[2], p[3], p[0]) - r[0]; if (td > d) d = td;
    }
    M->maxclear[t] = d;
}

static int is_solvent_accessible(const mole_complex *M, int t, double radius)
{
    double dx = M->center[3*t]   - M->vcenter[3*t];
    double dy = M->center[3*t+1] - M->vcenter[3*t+1];
    double dz = M->center[3*t+2] - M->vcenter[3*t+2];
    return M->maxclear[t] > 2*radius
        || (dx*dx + dy*dy + dz*dz) > radius*radius;
}

static int is_interior(const mole_complex *M, int t, double thr)
{
    return M->maxclear[t] < 2*thr;
}

/* Edge.Update. The clearance is measured only against the THREE atoms of the
   facet the two tetrahedra share - a local quantity, not a global minimum over
   all atoms. */
void mole_edge_update(mole_complex *M, int t, int k)
{
    int nb = M->tn[4*t+k], i;
    double cl = 1e300, len;
    if (nb < 0) return;
    for (i = 0; i < 4; i++) {
        int a;
        double d;
        if (i == k) continue;            /* the vertex opposite the shared face */
        a = M->tv[4*t+i];
        d = line_dist(&M->vcenter[3*t], &M->vcenter[3*nb], &M->axyz[3*a]) - M->arad[a];
        if (d < cl) cl = d;
    }
    len = dist3(&M->vcenter[3*t], &M->vcenter[3*nb]);
    M->eclear[4*t+k] = cl < 0 ? 0.0 : cl;
    M->elen[4*t+k]   = len;
    M->eweight[4*t+k] = cl < 0 ? 100000000.0 : len / (cl*cl + 0.000001);

    /* MOLE's DEFAULT weight function is VoronoiScale, not LengthAndRadius:
       TunnelWeightFunction.VoronoiScale is 0 and ComplexParameters sets it
       explicitly. It multiplies the ordinary weight by the squared
       circumcentre-to-centroid offset of BOTH tetrahedra, penalising routes
       through degenerate (sliver) tetrahedra.

       Getting this wrong is not subtle: with the plain weight our Dijkstra found
       a genuinely cheaper route than MOLE's own (0.3912 against 0.4788 on an
       identical graph), which is what exposed it. */
    {
        double sa = 0.0, sb = 0.0, d;
        int q;
        for (q = 0; q < 3; q++) {
            d = M->vcenter[3*t+q]  - M->center[3*t+q];  sa += d*d;
            d = M->vcenter[3*nb+q] - M->center[3*nb+q]; sb += d*d;
        }
        M->evweight[4*t+k] = M->eweight[4*t+k] * (sa + sb);
    }
}

/* Complex.ComputeDepth / ComputeDepthLength: distance from the nearest
   BOUNDARY vertex - hop count, or summed edge length. Run over the graph as it
   stands, so the order relative to vertex removal matters.

   MULTI-SOURCE. This used to run a separate relaxation from every boundary
   vertex into shared arrays, which computes the same minimum the slow way -
   O(boundary x V x E) - and it carried a fixed nt queue whose `tail` counted
   TOTAL enqueues rather than live occupancy. On hitting that cap it executed
   `head = 0; tail = 0`, DISCARDING every pending relaxation and abandoning
   that source, leaving vertices at their 1e300 / INT_MAX sentinel to be read
   as real depths. Simply unbounding the queue fixes the correctness bug and
   makes the cost unusable (measured: it hangs), because the per-source loop is
   then Bellman-Ford once per source.

   Seeding every boundary vertex at distance 0 removes both problems:
     hop count  - plain BFS. FIFO order makes depths non-decreasing, so each
                  vertex improves at most once and nt slots always suffice.
     path length- Dijkstra. Edge lengths are non-negative, so a popped vertex
                  is final and each directed edge yields at most one push;
                  5*nt + 8 slots therefore cannot overflow. */

typedef struct { double d; int v; } dh_item;

static void dh_push(dh_item *h, int *n, double d, int v)
{
    int i = (*n)++;
    h[i].d = d; h[i].v = v;
    while (i > 0) {
        int p = (i - 1) / 2;
        dh_item t;
        if (h[p].d <= h[i].d) break;
        t = h[p]; h[p] = h[i]; h[i] = t;
        i = p;
    }
}

static dh_item dh_pop(dh_item *h, int *n)
{
    dh_item top = h[0];
    int i = 0;
    h[0] = h[--(*n)];
    for (;;) {
        int l = 2*i + 1, r = l + 1, m = i;
        dh_item t;
        if (l < *n && h[l].d < h[m].d) m = l;
        if (r < *n && h[r].d < h[m].d) m = r;
        if (m == i) break;
        t = h[m]; h[m] = h[i]; h[i] = t;
        i = m;
    }
    return top;
}

static void compute_depth(mole_complex *M, int use_len)
{
    int i, k;
    for (i = 0; i < M->nt; i++) {
        if (!M->alive[i]) continue;
        if (use_len) M->depthlen[i] = M->boundary[i] ? 0.0 : 1e300;
        else         M->depth[i]    = M->boundary[i] ? 0 : 0x7fffffff;
    }
    if (!use_len) {
        int head = 0, tail = 0;
        int *q = malloc((size_t)(M->nt > 0 ? M->nt : 1) * sizeof(int));
        if (!q) return;
        for (i = 0; i < M->nt; i++)
            if (M->alive[i] && M->boundary[i]) q[tail++] = i;
        while (head < tail) {
            int c = q[head++];
            for (k = 0; k < 4; k++) {
                int v = M->tn[4*c+k];
                if (v < 0 || !M->alive[v]) continue;
                /* c is always a REACHED vertex (a source, or relaxed before it
                   was queued), so depth[c] is finite and cannot overflow. */
                if (M->depth[c] + 1 < M->depth[v]) {
                    M->depth[v] = M->depth[c] + 1;
                    q[tail++] = v;
                }
            }
        }
        free(q);
        return;
    }
    {
        int cap = 5 * (M->nt > 0 ? M->nt : 1) + 8;
        int hn = 0;
        dh_item *h = malloc((size_t)cap * sizeof(dh_item));
        if (!h) return;
        for (i = 0; i < M->nt; i++)
            if (M->alive[i] && M->boundary[i]) dh_push(h, &hn, 0.0, i);
        while (hn > 0) {
            dh_item it = dh_pop(h, &hn);
            int c = it.v;
            if (it.d > M->depthlen[c]) continue;   /* stale heap entry */
            for (k = 0; k < 4; k++) {
                int v = M->tn[4*c+k];
                double nd;
                if (v < 0 || !M->alive[v]) continue;
                nd = M->depthlen[c] + M->elen[4*c+k];
                if (nd < M->depthlen[v]) {
                    M->depthlen[v] = nd;
                    dh_push(h, &hn, nd, v);
                }
            }
        }
        free(h);
    }
}

/* Complex.RemoveShallowVertices. Not the procedure either paper describes; this
   follows the code. A vertex is "safe" if it neighbours the MinDepth+1 shell,
   propagated inward; anything else at depth <= MinDepth whose neighbours are all
   at most as deep is a surface ridge and goes. */
static void remove_shallow(mole_complex *M, int min_depth)
{
    char *safe = calloc((size_t)M->nt, 1);
    int i, k, d;
    if (!safe) return;
    for (i = 0; i < M->nt; i++) {
        if (!M->alive[i] || M->depth[i] != min_depth + 1) continue;
        for (k = 0; k < 4; k++) {
            int v = M->tn[4*i+k];
            if (v >= 0 && M->alive[v]) safe[v] = 1;
        }
    }
    /* Their loop takes a SNAPSHOT of the safe vertices at this depth before
       expanding, so anything that becomes safe during the pass is not itself
       expanded from until a later depth. Expanding as we go marks strictly more
       vertices safe and leaves 104 tetrahedra alive that MOLE removes. */
    {
        int *snap = xa_malloc((size_t)M->nt * sizeof(int));
        if (snap)
            for (d = min_depth; d >= 0; d--) {
                int ns = 0;
                for (i = 0; i < M->nt; i++)
                    if (M->alive[i] && safe[i] && M->depth[i] == d) snap[ns++] = i;
                for (i = 0; i < ns; i++)
                    for (k = 0; k < 4; k++) {
                        int v = M->tn[4*snap[i]+k];
                        if (v >= 0 && M->alive[v]) safe[v] = 1;
                    }
            }
        free(snap);
    }
    for (d = min_depth; d >= 0; d--) {
        int *kill = malloc((size_t)M->nt * sizeof(int));
        int nk = 0;
        if (!kill) break;
        for (i = 0; i < M->nt; i++) {
            int all_le = 1;
            if (!M->alive[i] || safe[i] || M->depth[i] != d) continue;
            for (k = 0; k < 4; k++) {
                int v = M->tn[4*i+k];
                if (v >= 0 && M->alive[v] && M->depth[v] > d) { all_le = 0; break; }
            }
            if (all_le) kill[nk++] = i;
        }
        for (i = 0; i < nk; i++) M->alive[kill[i]] = 0;
        free(kill);
    }
    free(safe);
}

/* Complex.RemoveOutsideAndInterior: peel inward from the convex hull while the
   probe still fits, then drop what is too narrow to hold a tunnel.
   IsBoundary is set on EVERY neighbour of a removed tetrahedron, whether or not
   that neighbour is itself removed - which is what makes the boundary layer one
   step wider than the removed set. */
static void remove_outside_and_interior(mole_complex *M, double probe, double interior)
{
    int *snap = xa_malloc((size_t)M->nt * sizeof(int));
    int *layer = malloc((size_t)M->nt * sizeof(int));
    char *seen = calloc((size_t)M->nt, 1);
    int ns = 0, i, k;
    if (!snap || !layer || !seen) return;

    for (i = 0; i < M->nt; i++)
        if (M->boundary[i]) { seen[i] = 1; M->alive[i] = 0; snap[ns++] = i; }

    while (ns > 0) {
        int nl = 0;
        for (i = 0; i < ns; i++)
            for (k = 0; k < 4; k++) {
                int of = M->tn[4*snap[i]+k];
                if (of < 0) continue;
                if (is_solvent_accessible(M, of, probe) && !seen[of]) {
                    seen[of] = 1; M->alive[of] = 0; layer[nl++] = of;
                }
                M->boundary[of] = 1;
            }
        memcpy(snap, layer, (size_t)nl * sizeof(int));
        ns = nl;
    }

    /* SurfaceCavity is snapshotted here, between the peel and the interior
       removal - it is not the same graph as the cavities. */
    M->n_surface = 0;
    for (i = 0; i < M->nt; i++) {
        if (M->surface) M->surface[i] = M->alive[i];
        if (M->alive[i]) M->n_surface++;
    }

    for (i = 0; i < M->nt; i++)
        if (M->alive[i] && is_interior(M, i, interior)) M->alive[i] = 0;

    free(snap); free(layer); free(seen);
}

/* Everything after the triangulation. Split out so MOLE's own tetrahedra can be
   substituted for ours - which is how the port is checked independently of the
   0.03% triangulation difference. */

/* cavityGraph.Vertices order, which mole_openings filters.
 *
 * AddVerticesAndEdgeRange(Triangulation.Edges) adds Source then Target, so a
 * tetrahedron takes its place at the first edge touching it. Edges come from
 * VoronoiMesh3D.Create walking cells in order, four neighbours each; a .NET
 * HashSet enumerates its slot array, so that walk is the edge order. Each
 * undirected edge is added once, at the smaller cell index - VoronoiEdge3D's
 * Equals is symmetric and rejects the duplicate. */
static void compute_vorder(mole_complex *M)
{
    int c, k, next = 0;
    for (c = 0; c < M->nt; c++) M->vorder[c] = -1;
    for (c = 0; c < M->nt; c++)
        for (k = 0; k < 4; k++) {
            int nb = M->tn[4*c+k];
            if (nb < 0 || nb < c) continue;   /* added already, at nb */
            if (M->vorder[c] < 0)  M->vorder[c]  = next++;
            if (M->vorder[nb] < 0) M->vorder[nb] = next++;
        }
    /* A cell with no finite neighbour never becomes a graph vertex; it cannot
       reach the opening code either, but leave it ordered last rather than -1
       so a comparison never sees a negative key. */
    for (c = 0; c < M->nt; c++) if (M->vorder[c] < 0) M->vorder[c] = next++;
}

int mole_build_from_tetra(mole_complex *M, const int *tv, const int *tn, int nt,
                          const double *axyz, const double *arad,
                          const mole_params *P)
{
    int i, k;
    memset(M, 0, sizeof(*M));
    M->nt = nt;
    M->axyz = axyz; M->arad = arad;
    M->tv = xa_malloc((size_t)nt*4*sizeof(int));
    M->tn = xa_malloc((size_t)nt*4*sizeof(int));
    M->center = xa_malloc((size_t)nt*3*sizeof(double));
    M->vcenter = xa_malloc((size_t)nt*3*sizeof(double));
    M->volume = xa_malloc((size_t)nt*sizeof(double));
    M->maxclear = xa_malloc((size_t)nt*sizeof(double));
    /* calloc, not malloc: ComputeDepth initialises only LIVE tetrahedra, so a
       dead one's depth is whatever the allocator handed back. Every read is
       guarded by alive[], so nothing used it - but it is still uninitialised
       memory, and it made the arrays non-deterministic, which the Tcl port's
       quantity-by-quantity comparison surfaced immediately. Zeroing changes no
       guarded read: verified byte-identical tunnel output. */
    M->depth = xa_calloc((size_t)nt, sizeof(int));
    M->depthlen = xa_calloc((size_t)nt, sizeof(double));
    M->eclear = xa_malloc((size_t)nt*4*sizeof(double));
    M->elen = xa_malloc((size_t)nt*4*sizeof(double));
    M->eweight = xa_malloc((size_t)nt*4*sizeof(double));
    M->evweight = xa_malloc((size_t)nt*4*sizeof(double));
    M->boundary = xa_calloc((size_t)nt, 1);
    M->alive = xa_malloc((size_t)nt);
    M->surface = xa_calloc((size_t)nt, 1);
    M->comp = xa_malloc((size_t)nt*sizeof(int));
    M->vorder = xa_malloc((size_t)nt*sizeof(int));
    if (!M->tv || !M->tn || !M->center || !M->vcenter || !M->volume || !M->maxclear
        || !M->depth || !M->depthlen || !M->eclear || !M->elen || !M->eweight || !M->evweight
        || !M->vorder || !M->surface
        || !M->boundary || !M->alive || !M->comp) return -1;

    memcpy(M->tv, tv, (size_t)nt*4*sizeof(int));
    memcpy(M->tn, tn, (size_t)nt*4*sizeof(int));
    for (i = 0; i < nt; i++) M->alive[i] = 1;
    for (i = 0; i < nt; i++)
        for (k = 0; k < 4; k++)
            if (M->tn[4*i+k] < 0) { M->boundary[i] = 1; break; }

    /* eclear is malloc'd, so the strict branch would read indeterminate bytes.
       Zero is the FAITHFUL value, not a safety default: it is what C# leaves
       Edge.Clearance at until the first live Edge.Update below. Do not "fix"
       this by calling mole_edge_update first - that is the line MOLE has
       commented out. */
    compute_vorder(M);
    if (P->strict_interior) memset(M->eclear, 0, (size_t)nt*4*sizeof(double));
    for (i = 0; i < nt; i++) tet_init(M, i, P->strict_interior);
    remove_outside_and_interior(M, P->probe_radius, P->interior_threshold);
    compute_depth(M, 0);
    for (i = 0; i < nt; i++) for (k = 0; k < 4; k++) mole_edge_update(M, i, k);
    remove_shallow(M, P->min_depth);
    compute_depth(M, 0);
    compute_depth(M, 1);
    return 0;
}

int mole_build(mole_complex *M, const dt_mesh *m, const double *axyz,
               const double *arad, const mole_params *P)
{
    int i, k, nt = 0;
    int *map = malloc((size_t)m->nt * sizeof(int));
    if (!map) return -1;
    for (i = 0; i < m->nt; i++)
        map[i] = (!m->t[i].dead && dt_is_finite(m, i)) ? nt++ : -1;

    memset(M, 0, sizeof(*M));
    M->nt = nt;
    M->axyz = axyz; M->arad = arad;
    M->tv = xa_malloc((size_t)nt*4*sizeof(int));
    M->tn = xa_malloc((size_t)nt*4*sizeof(int));
    M->center = xa_malloc((size_t)nt*3*sizeof(double));
    M->vcenter = xa_malloc((size_t)nt*3*sizeof(double));
    M->volume = xa_malloc((size_t)nt*sizeof(double));
    M->maxclear = xa_malloc((size_t)nt*sizeof(double));
    /* calloc, not malloc: ComputeDepth initialises only LIVE tetrahedra, so a
       dead one's depth is whatever the allocator handed back. Every read is
       guarded by alive[], so nothing used it - but it is still uninitialised
       memory, and it made the arrays non-deterministic, which the Tcl port's
       quantity-by-quantity comparison surfaced immediately. Zeroing changes no
       guarded read: verified byte-identical tunnel output. */
    M->depth = xa_calloc((size_t)nt, sizeof(int));
    M->depthlen = xa_calloc((size_t)nt, sizeof(double));
    M->eclear = xa_malloc((size_t)nt*4*sizeof(double));
    M->elen = xa_malloc((size_t)nt*4*sizeof(double));
    M->eweight = xa_malloc((size_t)nt*4*sizeof(double));
    M->evweight = xa_malloc((size_t)nt*4*sizeof(double));
    M->boundary = xa_calloc((size_t)nt, 1);
    M->alive = xa_malloc((size_t)nt);
    M->surface = xa_calloc((size_t)nt, 1);
    M->comp = xa_malloc((size_t)nt*sizeof(int));
    M->vorder = xa_malloc((size_t)nt*sizeof(int));
    if (!M->tv || !M->tn || !M->center || !M->vcenter || !M->volume || !M->maxclear
        || !M->depth || !M->depthlen || !M->eclear || !M->elen || !M->eweight || !M->evweight
        || !M->vorder || !M->surface
        || !M->boundary || !M->alive || !M->comp) { free(map); return -1; }

    for (i = 0; i < m->nt; i++) {
        int t = map[i];
        if (t < 0) continue;
        for (k = 0; k < 4; k++) {
            int nb = m->t[i].nb[k];
            M->tv[4*t+k] = m->t[i].v[k];
            M->tn[4*t+k] = (nb >= 0 && !m->t[nb].dead && dt_is_finite(m, nb)) ? map[nb] : -1;
        }
        M->alive[t] = 1;
    }
    free(map);

    /* A tetrahedron is boundary when it has fewer than four finite neighbours,
       i.e. it touches the infinite vertex. */
    for (i = 0; i < nt; i++)
        for (k = 0; k < 4; k++)
            if (M->tn[4*i+k] < 0) { M->boundary[i] = 1; break; }

    compute_vorder(M);
    if (P->strict_interior) memset(M->eclear, 0, (size_t)nt*4*sizeof(double));
    for (i = 0; i < nt; i++) tet_init(M, i, P->strict_interior);

    remove_outside_and_interior(M, P->probe_radius, P->interior_threshold);
    compute_depth(M, 0);
    for (i = 0; i < nt; i++) for (k = 0; k < 4; k++) mole_edge_update(M, i, k);
    remove_shallow(M, P->min_depth);
    compute_depth(M, 0);
    compute_depth(M, 1);
    return 0;
}

/* Connected components of the surviving graph, then the channel/void split.
   Channels touch the boundary and must clear both depth filters; voids do not
   touch it and must have positive volume and more than 20 tetrahedra. Both are
   ordered by descending volume, which is what fixes the cavity numbering. */
int mole_cavities(mole_complex *M, const mole_params *P,
                  mole_cavity **out, int *n_channel, int *n_void)
{
    int i, k, nc = 0, *stack = malloc((size_t)M->nt * sizeof(int));
    mole_cavity *cav = NULL;
    if (!stack) return -1;
    for (i = 0; i < M->nt; i++) M->comp[i] = -1;

    for (i = 0; i < M->nt; i++) {
        int sp = 0;
        if (!M->alive[i] || M->comp[i] >= 0) continue;
        cav = xa_realloc(cav, (size_t)(nc+1) * sizeof(*cav));
        memset(&cav[nc], 0, sizeof(cav[nc]));
        cav[nc].depth = 0; cav[nc].depth_length = 0.0;
        stack[sp++] = i; M->comp[i] = nc;
        while (sp > 0) {
            int c = stack[--sp];
            cav[nc].count++;
            cav[nc].volume += M->volume[c];
            if (M->boundary[c]) cav[nc].has_boundary = 1;
            if (M->depth[c] != 0x7fffffff && M->depth[c] > cav[nc].depth)
                cav[nc].depth = M->depth[c];
            if (M->depthlen[c] < 1e299 && M->depthlen[c] > cav[nc].depth_length)
                cav[nc].depth_length = M->depthlen[c];
            for (k = 0; k < 4; k++) {
                int v = M->tn[4*c+k];
                if (v >= 0 && M->alive[v] && M->comp[v] < 0) { M->comp[v] = nc; stack[sp++] = v; }
            }
        }
        nc++;
    }
    free(stack);

    *n_channel = 0; *n_void = 0;
    for (i = 0; i < nc; i++) {
        if (cav[i].has_boundary) {
            if (cav[i].depth_length > P->min_depth_length && cav[i].depth > P->min_depth)
                (*n_channel)++;
        } else if (cav[i].volume > 0 && cav[i].count > 20) {
            (*n_void)++;
        }
    }
    *out = cav;
    return nc;
}

void mole_free(mole_complex *M)
{
    free(M->tv); free(M->tn); free(M->center); free(M->vcenter);
    free(M->volume); free(M->maxclear); free(M->depth); free(M->depthlen);
    free(M->eclear); free(M->elen); free(M->eweight); free(M->evweight);
    free(M->vorder);
    free(M->boundary); free(M->alive); free(M->surface); free(M->comp);
    memset(M, 0, sizeof(*M));
}
