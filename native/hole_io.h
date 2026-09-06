/* hole_io.h - the one reader for HOLE's .sph and for PDB ATOM records.
 *
 * Every native tool used to slice the fixed PDB columns itself, and the
 * parsers drifted: one knew about HOLE's LAST-REC-END cutter marker, the
 * next did not. This iterator hands each ATOM/HETATM record out with every
 * column already parsed and the marker resolved; the caller keeps only its
 * own filter.
 *
 *   hio_reader r; hio_rec a;
 *   if (!hio_open(&r, path)) ...;
 *   while (hio_next(&r, &a)) { ... }
 *   hio_close(&r);
 *
 * Numeric columns parse the way Tcl's [string range] + [string is double]
 * do: blank-trimmed, whole field or nothing (the *_ok flags). Short lines
 * read as blank. a.line is the raw record with CR/LF stripped; it stays
 * valid until the next hio_next. */
#ifndef HOLE_IO_H
#define HOLE_IO_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    double x, y, z;      /* cols 30-37, 38-45, 46-53 */
    double occ, beta;    /* cols 54-59, 60-65 */
    int resseq;          /* cols 22-25 */
    int xyz_ok, occ_ok, beta_ok, resseq_ok;
    int marked;          /* a LAST-REC-END line followed this record */
    char name[5];        /* cols 12-15, raw */
    char resn[4];        /* cols 17-19, raw */
    char chain, icode;   /* col 21, col 26 */
    char segid[5];       /* cols 72-75, raw, blank when short */
    const char *line;
    int len;
} hio_rec;

typedef struct {
    FILE *f;
    char buf[2][1024];
    int cur, have;
    hio_rec pend;
} hio_reader;

static void hio_field(const char *line, int len, int a, int b, char *out, int outsz) {
    int n = b - a + 1, k = 0;
    for (int i = 0; i < n && k < outsz - 1; i++) {
        int p = a + i;
        out[k++] = (p >= 0 && p < len) ? line[p] : ' ';
    }
    out[k] = 0;
    while (k > 0 && (out[k-1] == ' ' || out[k-1] == '\t')) out[--k] = 0;
    int s = 0; while (out[s] == ' ') s++;
    if (s) memmove(out, out + s, strlen(out + s) + 1);
}

static int hio_num(const char *line, int len, int a, int b, double *v) {
    char t[32]; hio_field(line, len, a, b, t, sizeof t);
    *v = 0.0;
    if (!*t) return 0;
    char *end; *v = strtod(t, &end);
    while (*end == ' ') end++;
    return *end == 0;
}

static void hio_raw(const char *line, int len, int a, int n, char *out) {
    for (int i = 0; i < n; i++) out[i] = (a + i < len) ? line[a + i] : ' ';
    out[n] = 0;
}

static int hio_is_atom(const char *line) {
    return !strncmp(line, "ATOM  ", 6) || !strncmp(line, "HETATM", 6);
}

static void hio_parse(hio_rec *r, char *line) {
    int len = (int)strlen(line);
    while (len && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = 0;
    memset(r, 0, sizeof *r);
    r->line = line; r->len = len;
    int ox = hio_num(line, len, 30, 37, &r->x), oy = hio_num(line, len, 38, 45, &r->y),
        oz = hio_num(line, len, 46, 53, &r->z);
    r->xyz_ok = ox && oy && oz;
    r->occ_ok = hio_num(line, len, 54, 59, &r->occ);
    r->beta_ok = hio_num(line, len, 60, 65, &r->beta);
    double rs; r->resseq_ok = hio_num(line, len, 22, 25, &rs); r->resseq = (int)rs;
    hio_raw(line, len, 12, 4, r->name);
    hio_raw(line, len, 17, 3, r->resn);
    r->chain = len > 21 ? line[21] : ' ';
    r->icode = len > 26 ? line[26] : ' ';
    hio_raw(line, len, 72, 4, r->segid);
}

static int hio_open(hio_reader *r, const char *path) {
    memset(r, 0, sizeof *r);
    r->f = fopen(path, "r");
    return r->f != NULL;
}

/* The marker follows its record, so records are handed out one line late. */
static int hio_next(hio_reader *r, hio_rec *out) {
    for (;;) {
        char *nb = r->buf[1 - r->cur];
        if (!fgets(nb, sizeof r->buf[0], r->f)) {
            if (!r->have) return 0;
            *out = r->pend; r->have = 0; return 1;
        }
        if (!strncmp(nb, "LAST-REC-END", 12)) { if (r->have) r->pend.marked = 1; continue; }
        if (!hio_is_atom(nb)) continue;
        hio_rec nr; hio_parse(&nr, nb);
        r->cur = 1 - r->cur;
        if (r->have) { *out = r->pend; r->pend = nr; return 1; }
        r->pend = nr; r->have = 1;
    }
}

static void hio_close(hio_reader *r) { if (r->f) fclose(r->f); r->f = NULL; }
#endif
