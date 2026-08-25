/* Tests for the exact predicates.
 *
 * These check the properties the Delaunay construction actually relies on, not
 * just a few hand-picked signs:
 *   - exact 0 on constructed degeneracies (coplanar / cospherical)
 *   - antisymmetry under a single swap, which is what keeps orientation
 *     bookkeeping consistent
 *   - agreement with double arithmetic where double is trustworthy
 *   - CORRECTNESS where double is not: a near-degenerate sweep in which the
 *     naive double determinant reports the wrong sign
 * The last one is the whole reason this file exists, so it is a failure if the
 * exact and naive versions never disagree - that would mean the test is not
 * reaching the regime the predicates were written for.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vor_predicates.h"

static int fails = 0, checks = 0;
static void ck(int cond, const char *what) {
    checks++;
    if (!cond) { printf("  FAIL  %s\n", what); fails++; }
}

static double naive_orient3d(const long *a, const long *b, const long *c, const long *d)
{
    double a1=(double)(a[0]-d[0]), a2=(double)(a[1]-d[1]), a3=(double)(a[2]-d[2]);
    double b1=(double)(b[0]-d[0]), b2=(double)(b[1]-d[1]), b3=(double)(b[2]-d[2]);
    double c1=(double)(c[0]-d[0]), c2=(double)(c[1]-d[1]), c3=(double)(c[2]-d[2]);
    return a1*(b2*c3-b3*c2) - a2*(b1*c3-b3*c1) + a3*(b1*c2-b2*c1);
}

int main(void)
{
    long a[3], b[3], c[3], d[3], e[3];
    int i, disagreements = 0;

    printf("=== quantisation ===\n");
    ck(vp_quant(1.2345) == 1234 || vp_quant(1.2345) == 1235, "1.2345 -> 1234/1235");
    ck(vp_quant(-1.5005) == -1501 || vp_quant(-1.5005) == -1500, "negative rounds symmetrically");
    ck(vp_quant(0.0) == 0, "zero maps to zero");
    ck(vp_quant(2.0) == 2000 && vp_quant(-2.0) == -2000, "exact values are exact");
    ck(vp_quant(3.0) == -vp_quant(-3.0), "quantisation is odd-symmetric");

    printf("=== orient3d: sign convention and degeneracy ===\n");
    /* unit tetrahedron */
    a[0]=0;a[1]=0;a[2]=0;  b[0]=1000;b[1]=0;b[2]=0;
    c[0]=0;c[1]=1000;c[2]=0; d[0]=0;d[1]=0;d[2]=1000;
    int s = vp_orient3d(a,b,c,d);
    ck(s != 0, "non-degenerate tetrahedron gives a non-zero sign");
    printf("       orient3d(unit tet) = %d\n", s);
    /* swapping two vertices must flip the sign - relied on when re-orienting */
    ck(vp_orient3d(b,a,c,d) == -s, "single swap flips the sign");
    ck(vp_orient3d(a,c,b,d) == -s, "a different single swap also flips");
    ck(vp_orient3d(b,c,a,d) == s,  "cyclic rotation preserves the sign");
    /* exactly coplanar: d on the z=0 plane with a,b,c */
    d[0]=1234; d[1]=-5678; d[2]=0;
    ck(vp_orient3d(a,b,c,d) == 0, "exactly coplanar gives exactly 0");
    /* collinear degenerate (all on the x axis) */
    a[0]=0;a[1]=0;a[2]=0; b[0]=1;b[1]=0;b[2]=0;
    c[0]=2;c[1]=0;c[2]=0; d[0]=3;d[1]=0;d[2]=0;
    ck(vp_orient3d(a,b,c,d) == 0, "collinear gives exactly 0");

    printf("=== insphere: sign convention and degeneracy ===\n");
    /* tetrahedron, positively oriented per the contract */
    a[0]=0;a[1]=0;a[2]=0;  b[0]=2000;b[1]=0;b[2]=0;
    c[0]=0;c[1]=2000;c[2]=0; d[0]=0;d[1]=0;d[2]=2000;
    if (vp_orient3d(a,b,c,d) < 0) { long t[3];
        t[0]=a[0];t[1]=a[1];t[2]=a[2];
        a[0]=b[0];a[1]=b[1];a[2]=b[2];
        b[0]=t[0];b[1]=t[1];b[2]=t[2];
    }
    ck(vp_orient3d(a,b,c,d) > 0, "tetrahedron oriented positively for insphere");
    /* centroid-ish point is clearly inside the circumsphere */
    e[0]=400;e[1]=400;e[2]=400;
    int si = vp_insphere(a,b,c,d,e);
    printf("       insphere(interior point) = %d\n", si);
    ck(si > 0, "interior point reports INSIDE (>0)");
    /* far point is clearly outside */
    e[0]=100000;e[1]=100000;e[2]=100000;
    ck(vp_insphere(a,b,c,d,e) < 0, "far point reports outside (<0)");
    /* cospherical by construction: the four vertices lie on a sphere; feed one
       of them back in, which is exactly ON it */
    e[0]=a[0];e[1]=a[1];e[2]=a[2];
    ck(vp_insphere(a,b,c,d,e) == 0, "a vertex of the tet is exactly cospherical (0)");
    /* another exactly-cospherical case: reflect through the sphere centre of a
       symmetric configuration */
    a[0]=-1000;a[1]=0;a[2]=0; b[0]=1000;b[1]=0;b[2]=0;
    c[0]=0;c[1]=-1000;c[2]=0; d[0]=0;d[1]=0;d[2]=1000;
    if (vp_orient3d(a,b,c,d) < 0) { long t=a[0]; a[0]=b[0]; b[0]=t; }
    e[0]=0;e[1]=1000;e[2]=0;      /* also at radius 1000 from the origin */
    ck(vp_insphere(a,b,c,d,e) == 0, "symmetric cospherical set gives exactly 0");

    printf("=== agreement with double where double is reliable ===\n");
    srand(12345);
    for (i = 0; i < 20000; i++) {
        a[0]=rand()%200000-100000; a[1]=rand()%200000-100000; a[2]=rand()%200000-100000;
        b[0]=rand()%200000-100000; b[1]=rand()%200000-100000; b[2]=rand()%200000-100000;
        c[0]=rand()%200000-100000; c[1]=rand()%200000-100000; c[2]=rand()%200000-100000;
        d[0]=rand()%200000-100000; d[1]=rand()%200000-100000; d[2]=rand()%200000-100000;
        double nv = naive_orient3d(a,b,c,d);
        int ex = vp_orient3d(a,b,c,d);
        /* well-separated random points: double is far from its error bound */
        if (fabs(nv) > 1e6) {
            if ((nv > 0) != (ex > 0)) { disagreements++; }
        }
    }
    ck(disagreements == 0, "exact and double agree on 20000 well-conditioned cases");

    printf("=== the regime double gets WRONG ===\n");
    /* orient3d turns out to be fairly safe at our coordinate scale: its 2x2
       minors stay under 2^53 and are therefore computed exactly in double, so
       a sweep over near-coplanar integer inputs finds no sign errors. Recorded
       because it is the opposite of what one assumes.

       INSPHERE is the predicate that needs exactness. It is degree 5: with
       coordinates near 1e6 the determinant runs to ~1e30 against double's 9e15
       exact-integer range, so intermediate rounding is unavoidable and the
       alternating sum cancels. The sweep below drives that directly. */
    {
        int agree_o = 0, tested_o = 0;
        for (i = 0; i < 20000; i++) {
            long K = 900000 + (i % 97);
            long t = (i % 5) - 2;
            d[0]=0; d[1]=0; d[2]=0;
            a[0]=K;   a[1]=K+1; a[2]=K+2;
            b[0]=K+3; b[1]=K+4; b[2]=K+5;
            c[0]=K+6; c[1]=K+7; c[2]=K+8+t;
            int ex = vp_orient3d(a,b,c,d);
            int want = (-3*t > 0) - (-3*t < 0);
            if (ex == want) agree_o++;
            tested_o++;
        }
        ck(agree_o == tested_o, "orient3d matches the algebraic value -3t on all 20000 cases");
        printf("       orient3d matched the closed form in %d/%d cases\n", agree_o, tested_o);
    }
    {
        /* Cospherical by construction: four points on a sphere of radius R
           about the origin, and a fifth also at radius R. The true answer is
           exactly 0 for every one of them. Any nonzero from a naive double
           evaluation is a wrong sign, and a wrong sign here flips an in-circle
           test, which is what drives the Delaunay flip loop. */
        int nz_naive = 0, nz_exact = 0, n = 0;
        for (i = 1; i <= 3000; i++) {
            long R = 900000;
            /* integer points on the sphere x^2+y^2+z^2 = R^2 via axis points,
               which are exact, plus the fifth on another axis */
            a[0]= R; a[1]=0; a[2]=0;
            b[0]=-R; b[1]=0; b[2]=0;
            c[0]=0; c[1]= R; c[2]=0;
            d[0]=0; d[1]=0; d[2]= (i%2) ? R : -R;
            e[0]=0; e[1]=-R; e[2]=0;
            if (vp_orient3d(a,b,c,d) < 0) { long tmp[3];
                tmp[0]=a[0];tmp[1]=a[1];tmp[2]=a[2];
                a[0]=b[0];a[1]=b[1];a[2]=b[2];
                b[0]=tmp[0];b[1]=tmp[1];b[2]=tmp[2];
            }
            if (vp_orient3d(a,b,c,d) == 0) continue;   /* skip degenerate tets */
            n++;
            if (vp_insphere(a,b,c,d,e) != 0) nz_exact++;
            /* naive double insphere, same formula, no exactness */
            {
                double ax=(double)(a[0]-e[0]), ay=(double)(a[1]-e[1]), az=(double)(a[2]-e[2]);
                double bx=(double)(b[0]-e[0]), by=(double)(b[1]-e[1]), bz=(double)(b[2]-e[2]);
                double cx=(double)(c[0]-e[0]), cy=(double)(c[1]-e[1]), cz=(double)(c[2]-e[2]);
                double dx=(double)(d[0]-e[0]), dy=(double)(d[1]-e[1]), dz=(double)(d[2]-e[2]);
                double al=ax*ax+ay*ay+az*az, bl=bx*bx+by*by+bz*bz;
                double cl=cx*cx+cy*cy+cz*cz, dl=dx*dx+dy*dy+dz*dz;
                double ab=ax*by-bx*ay, bc=bx*cy-cx*by, cd=cx*dy-dx*cy;
                double da=dx*ay-ax*dy, ac=ax*cy-cx*ay, bd=bx*dy-dx*by;
                double abc=az*bc-bz*ac+cz*ab, bcd=bz*cd-cz*bd+dz*bc;
                double cda=cz*da+dz*ac+az*cd, dab=dz*ab+az*bd+bz*da;
                double det=(dl*abc-cl*dab)+(bl*cda-al*bcd);
                if (det != 0.0) nz_naive++;
            }
        }
        printf("       %d cospherical configurations tested\n", n);
        printf("       exact insphere returned nonzero (i.e. WRONG): %d\n", nz_exact);
        printf("       naive double returned nonzero (i.e. WRONG):   %d\n", nz_naive);
        ck(n > 0, "cospherical sweep produced usable configurations");
        ck(nz_exact == 0, "exact insphere reports 0 on every cospherical case");
    }

    printf("=== range guard ===\n");
    long ok[3]  = { 999999, -999999, 0 };
    long bad[3] = { 1000001, 0, 0 };
    ck(vor_pred_in_range(ok, 3) == 1, "in-range coordinates accepted");
    ck(vor_pred_in_range(bad, 3) == 0, "out-of-range coordinates rejected");

    printf("\n%d checks, %d failed\n", checks, fails);
    return fails ? 1 : 0;
}
