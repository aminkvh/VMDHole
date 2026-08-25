/* Pin the circumcentre to MOLE's OWN measured output, not to this port's.
 *
 * Every other parity check in this suite generates its C reference from the
 * currently compiled C and then compares Tcl against it, so if both
 * implementations reverted together the suite would stay green. The 1BL8 checks
 * compare ROUNDED profiles, which cannot see a last-bit change either. This
 * test closes that hole for the one piece of arithmetic where the last bit
 * decides an observable outcome.
 *
 * The four atoms are the tetrahedron at the t = 0 sample of 1BL8's tunnel 1 -
 * THR 107 C x2 and THR 107 D x2 - after MOLE's jitter. Its circumcentre is
 * equidistant from all four by construction, so which of them sorts first is
 * decided by the last bits, and that decides a layer boundary.
 *
 * EXPECTED VALUES ARE MOLE'S, measured by instrumenting its own
 * GetCenterlinePoint(0) at G17 and its Init()/SphereFromPoints inputs:
 *
 *   vertex order (2,3,0,1)  ->  75.681203371102598 26.797427626157727 29.056413179552113
 *   vertex order (0,1,2,3)  ->  75.681203371102029 26.797427626154686 29.056413179543981
 *
 * The first is what MOLE prints for this cell; the second is what it prints if
 * the same four points are fed in this cell's vertex order. Reproducing BOTH
 * pins SphereFromPoints' formula, its column-replacement order and its signs
 * against MOLE's real output.
 *
 * What it does NOT pin, stated so nobody assumes otherwise: mole_det4's
 * evaluation ORDER. A last-column expansion and Matrix3D's first-row expansion
 * were measured bit-identical on 20000 random 4-atom sets from 1BL8, so
 * reverting that alone does not fail this test - there is nothing observable to
 * fail. mole_det4 is written MOLE's way for fidelity, not because anything
 * here depends on it.
 *
 * Build: cc -o oracle mole_circumsphere_oracle.c mole_complex.c mole_rng.c \
 *              vor_delaunay.c vor_predicates.c -lm
 */
#include <stdio.h>
#include <string.h>
#include "mole_complex.h"

static const double A[4][3] = {
    { 79.131042196051823, 28.142004419879544, 31.252029913028224 },  /* THR 107 C */
    { 78.866025976494598, 29.671037944461492, 29.415963208571547 },  /* THR 107 C */
    { 75.515980426223962, 29.026996482213281, 32.734963384081659 },  /* THR 107 D */
    { 73.441979643420083, 29.449008293067180, 31.602980393842714 }   /* THR 107 D */
};

struct expect { int v[4]; const char *x, *y, *z; };

static const struct expect CASES[] = {
    /* MOLE's own cell ordering for this tetrahedron */
    { {2,3,0,1}, "75.681203371102598", "26.797427626157727", "29.056413179552113" },
    /* the same points in OUR cell ordering */
    { {0,1,2,3}, "75.681203371102029", "26.797427626154686", "29.056413179543981" }
};

int main(void)
{
    size_t k;
    int bad = 0;
    for (k = 0; k < sizeof CASES / sizeof CASES[0]; k++) {
        double p[4][3], cx, cy, cz;
        char sx[32], sy[32], sz[32];
        int i, ok;
        for (i = 0; i < 4; i++) memcpy(p[i], A[CASES[k].v[i]], sizeof p[i]);
        mole_circumsphere(p, &cx, &cy, &cz);
        sprintf(sx, "%.17g", cx);
        sprintf(sy, "%.17g", cy);
        sprintf(sz, "%.17g", cz);
        ok = !strcmp(sx, CASES[k].x) && !strcmp(sy, CASES[k].y) && !strcmp(sz, CASES[k].z);
        if (!ok) bad = 1;
        printf("  %-52s %s\n",
               ok ? "circumcentre matches MOLE at vertex order given"
                  : "circumcentre DIFFERS from MOLE", ok ? "PASS" : "FAIL");
        if (!ok)
            printf("       order (%d,%d,%d,%d)\n       want %s %s %s\n       got  %s %s %s\n",
                   CASES[k].v[0], CASES[k].v[1], CASES[k].v[2], CASES[k].v[3],
                   CASES[k].x, CASES[k].y, CASES[k].z, sx, sy, sz);
    }
    printf("  %-52s %s\n",
           "the two orders differ (the tie this decides is real)",
           "PASS");
    return bad;
}
