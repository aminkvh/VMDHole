/* Regression test for the mole_auto_origins buffer overflow (found in code
 * review, item 3).
 *
 * The caller passes a fixed `int origins[MOLE_MAX_ORIGINS]`; handing
 * mole_auto_origins an UNCLAMPED user value (the GUI accepted any nonnegative
 * integer), and mole_auto_origins' first write was unconditional:
 *
 *     out[n++] = cand[0];        // before any n < max_origins test
 *
 * so it stored an element even for max_origins == 0. Two checks here:
 *
 *   CANARY  call with a huge max_origins and verify nothing is written past
 *           MOLE_MAX_ORIGINS. Needs a cavity with enough qualifying maxima to
 *           actually reach the cap, which ordinary fixtures do not have - so
 *           this one is a guard against future regressions rather than a
 *           reproduction of the original bug.
 *   ZERO    call with max_origins == 0 and verify NOTHING is written. This one
 *           IS deterministic on any fixture: it fails against the pre-fix
 *           routine on the very first cavity, no 17-maxima structure needed.
 *
 * Usage: mole_maxorigins_test ATOMS.txt
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mole_rng.h"
#include "mole_complex.h"
#include "mole_tunnel.h"

#define CANARY 0x5AFEBEEF
#define SLACK  32

static mole_atoms A;
static double piv[3*MOLE_MAXATOM], rad[MOLE_MAXATOM];
static int    pivsrc[MOLE_MAXATOM];

int main(int argc, char **argv)
{
    dt_mesh m; mole_complex M;
    mole_params P = {0};
    mole_cavity *cav = NULL;
    int np, i, c, nc, nch, nvd;
    int fail = 0, tested = 0;
    int buf[MOLE_MAX_ORIGINS + SLACK];

    if (argc < 2) { fprintf(stderr, "usage: mole_maxorigins_test ATOMS.txt\n"); return 2; }
    if (mole_read_atoms(argv[1], &A, MOLE_MAXATOM) < 4) {
        fprintf(stderr, "cannot read %s\n", argv[1]); return 2;
    }
    np = mole_pivots(&A, piv, rad, pivsrc);
    if (dt_build(&m, piv, np) != 0) { fprintf(stderr, "dt_build failed\n"); return 1; }
    P.probe_radius = 3.0; P.interior_threshold = 1.25;
    P.min_depth = 8; P.min_depth_length = 5.0;
    if (mole_build(&M, &m, piv, rad, &P) != 0) { fprintf(stderr, "mole_build failed\n"); return 1; }
    nc = mole_cavities(&M, &P, &cav, &nch, &nvd);

    for (c = 0; c < nc; c++) {
        int n;
        if (!(cav[c].has_boundary && cav[c].depth_length > P.min_depth_length
              && cav[c].depth > P.min_depth)) continue;
        tested++;

        /* --- CANARY: an absurd max_origins must not write past the array --- */
        for (i = 0; i < MOLE_MAX_ORIGINS + SLACK; i++) buf[i] = CANARY;
        n = mole_auto_origins(&M, c, 10.0, 1000000, buf);
        if (n > MOLE_MAX_ORIGINS) {
            printf("  FAIL  cavity %d: returned %d origins, cap is %d\n",
                   c, n, MOLE_MAX_ORIGINS);
            fail++;
        }
        for (i = MOLE_MAX_ORIGINS; i < MOLE_MAX_ORIGINS + SLACK; i++)
            if (buf[i] != CANARY) {
                printf("  FAIL  cavity %d: wrote past the array at index %d\n", c, i);
                fail++;
                break;
            }

        /* --- ZERO: max_origins == 0 must write nothing at all --- */
        for (i = 0; i < MOLE_MAX_ORIGINS + SLACK; i++) buf[i] = CANARY;
        n = mole_auto_origins(&M, c, 10.0, 0, buf);
        if (n != 0) {
            printf("  FAIL  cavity %d: max_origins=0 returned %d\n", c, n);
            fail++;
        }
        if (buf[0] != CANARY) {
            printf("  FAIL  cavity %d: max_origins=0 still wrote out[0]\n", c);
            fail++;
        }
    }

    if (tested == 0) {
        printf("  FAIL  no qualifying cavity in %s - this test checked nothing\n", argv[1]);
        return 1;
    }
    if (fail == 0)
        printf("  PASS  mole_auto_origins respects its buffer on %d cavity(ies)\n", tested);
    return fail ? 1 : 0;
}
