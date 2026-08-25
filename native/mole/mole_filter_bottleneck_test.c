/* mole_filter_bottleneck's n<1 branch: MOLE keeps a fully-degenerate tunnel
 * UNCONDITIONALLY, derived from FilterBottleneck's own NaN comparisons - see
 * the comment on the branch in mole_tunnel.c. This exercises exactly that
 * branch directly, since no real structure reaches it (shortest candidate
 * profile across the fixtures is 3.47 A; see check 8b in
 * test_mole_tcl_port.sh for that measurement).
 *
 * A degenerate profile has length < 1/density. Any length under that with
 * garbage-but-finite spline data must return 1 (keep), regardless of
 * bottleneck, tolerance, or what the spline would evaluate to - because MOLE
 * never gets far enough to look at it. */
#include <stdio.h>
#include <string.h>
#include "mole_tunnel.h"

static void make_profile(mole_tunnel_profile *p, double length)
{
    static double t[2] = { 0.0, 1.0 };
    /* A radius BELOW the bottleneck and coordinates that would accumulate
       real distance if evaluated - if the branch under test ever regresses
       to actually running the loop, these values would make it reject. */
    static double lowr[2] = { 0.01, 0.01 };
    static double xyz[2] = { 0.0, 100.0 };
    memset(p, 0, sizeof *p);
    p->length = length;
    mole_spline_init(&p->sr, t, lowr, 2);
    mole_spline_init(&p->sx, t, xyz, 2);
    mole_spline_init(&p->sy, t, xyz, 2);
    mole_spline_init(&p->sz, t, xyz, 2);
}

int main(void)
{
    mole_tunnel_profile p;
    int bad = 0, i;
    /* length * density < 1 at every density any call site uses (2, 6, 8). */
    double lens[] = { 0.01, 0.05, 0.1, 0.12 };
    for (i = 0; i < (int)(sizeof lens / sizeof lens[0]); i++) {
        int r;
        make_profile(&p, lens[i]);
        /* bottleneck 5.0 (above the 0.01 radius) and tolerance BOTH 0 and
           nonzero - FilterBottleneck's two branches, both derived to
           unconditionally keep. */
        r = mole_filter_bottleneck(&p, 5.0, 0.0, 8.0);
        printf("  length=%.2f tolerance=0.0  -> %s\n", lens[i], r ? "keep" : "REJECT");
        if (!r) bad++;
        r = mole_filter_bottleneck(&p, 5.0, 2.0, 8.0);
        printf("  length=%.2f tolerance=2.0  -> %s\n", lens[i], r ? "keep" : "REJECT");
        if (!r) bad++;
        mole_profile_free(&p);
    }
    if (bad) {
        printf("  FAIL: %d of %d degenerate cases were rejected - MOLE keeps all of them\n",
               bad, (int)(sizeof lens / sizeof lens[0]) * 2);
        return 1;
    }
    printf("  PASS: every degenerate profile kept, matching MOLE's derived verdict\n");
    return 0;
}
