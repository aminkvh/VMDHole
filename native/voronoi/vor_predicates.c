/* Exact orient3d / insphere on quantised integer coordinates.
   See vor_predicates.h for why this is integer arithmetic and not floating
   point, and for the range argument. */
#include <math.h>
#include "vor_predicates.h"

long vp_quant(double x)
{
    double s = x * VP_SCALE;
    /* Round half away from zero, so the grid is symmetric about the origin and
       does not depend on the FP rounding mode. floor(x+0.5) would bias one way
       for negatives and make the predicates asymmetric under reflection. */
    return (long)(s >= 0.0 ? floor(s + 0.5) : ceil(s - 0.5));
}

int vor_pred_in_range(const long *p, int n3)
{
    int i;
    for (i = 0; i < n3; i++) {
        if (p[i] > VP_MAX_COORD || p[i] < -VP_MAX_COORD) return 0;
    }
    return 1;
}

static int vp_sign(vp_int v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

/* 256-bit signed accumulation, for insphere's final combination only.
 *
 * insphere's minors and lifts each fit vp_int comfortably; it is their PRODUCT,
 * the degree-5 term, that can overflow. At the default grid nothing does. At the
 * finer grid the MOLE port needs, the super-tetrahedron sits at 20x the data
 * extent and the product reaches ~9e43 against __int128's 1.7e38, which
 * silently returns a wrong SIGN - and a wrong sign gives corrupt topology, not
 * a slightly wrong number. dt_build simply failed.
 *
 * Widening cannot change any sign that was already being computed correctly, so
 * this is unconditional rather than a second code path to keep in step. */
typedef struct { unsigned __int128 lo; __int128 hi; } vp_i256;

static vp_i256 vp_mul_wide(vp_int a, vp_int b)
{
    /* Magnitudes multiplied as four 64x64 limb products, sign applied after. */
    int neg = (a < 0) ^ (b < 0);
    unsigned __int128 ua = (unsigned __int128)(a < 0 ? -a : a);
    unsigned __int128 ub = (unsigned __int128)(b < 0 ? -b : b);
    unsigned long long a0 = (unsigned long long)ua, a1 = (unsigned long long)(ua >> 64);
    unsigned long long b0 = (unsigned long long)ub, b1 = (unsigned long long)(ub >> 64);
    unsigned __int128 p00 = (unsigned __int128)a0 * b0;
    unsigned __int128 p01 = (unsigned __int128)a0 * b1;
    unsigned __int128 p10 = (unsigned __int128)a1 * b0;
    unsigned __int128 p11 = (unsigned __int128)a1 * b1;
    unsigned __int128 mid = (p00 >> 64) + (unsigned long long)p01 + (unsigned long long)p10;
    vp_i256 r;
    r.lo = (p00 & 0xFFFFFFFFFFFFFFFFULL) | (mid << 64);
    r.hi = (__int128)(p11 + (p01 >> 64) + (p10 >> 64) + (mid >> 64));
    if (neg) {   /* two's complement negate across both limbs */
        r.lo = ~r.lo + 1;
        r.hi = ~r.hi + (r.lo == 0 ? 1 : 0);
    }
    return r;
}

static vp_i256 vp_add256(vp_i256 x, vp_i256 y)
{
    vp_i256 r;
    r.lo = x.lo + y.lo;
    r.hi = x.hi + y.hi + (r.lo < x.lo ? 1 : 0);
    return r;
}

static int vp_sign256(vp_i256 v)
{
    if (v.hi != 0) return v.hi > 0 ? 1 : -1;
    return v.lo != 0 ? 1 : 0;
}

int vp_orient3d(const long *a, const long *b, const long *c, const long *d)
{
    /* |ax-dx ay-dy az-dz|
       |bx-dx by-dy bz-dz|
       |cx-dx cy-dy cz-dz|
       Differences first, so every term stays within the range argued in the
       header; expanding about the origin instead would square the magnitudes. */
    vp_int a1 = (vp_int)(a[0]-d[0]), a2 = (vp_int)(a[1]-d[1]), a3 = (vp_int)(a[2]-d[2]);
    vp_int b1 = (vp_int)(b[0]-d[0]), b2 = (vp_int)(b[1]-d[1]), b3 = (vp_int)(b[2]-d[2]);
    vp_int c1 = (vp_int)(c[0]-d[0]), c2 = (vp_int)(c[1]-d[1]), c3 = (vp_int)(c[2]-d[2]);
    return vp_sign( a1*(b2*c3 - b3*c2)
                  - a2*(b1*c3 - b3*c1)
                  + a3*(b1*c2 - b2*c1) );
}

int vp_insphere(const long *a, const long *b, const long *c,
                const long *d, const long *e)
{
    /* |ax-ex ay-ey az-ez |a-e|^2|
       |bx-ex ...                |
       |cx-ex ...                |
       |dx-ex ...                |
       Lifting onto the paraboloid: the 4x4 determinant is positive when e is
       inside the circumsphere of a,b,c,d, PROVIDED abcd is positively
       oriented - which is the caller's contract, stated in the header. */
    vp_int ax = (vp_int)(a[0]-e[0]), ay = (vp_int)(a[1]-e[1]), az = (vp_int)(a[2]-e[2]);
    vp_int bx = (vp_int)(b[0]-e[0]), by = (vp_int)(b[1]-e[1]), bz = (vp_int)(b[2]-e[2]);
    vp_int cx = (vp_int)(c[0]-e[0]), cy = (vp_int)(c[1]-e[1]), cz = (vp_int)(c[2]-e[2]);
    vp_int dx = (vp_int)(d[0]-e[0]), dy = (vp_int)(d[1]-e[1]), dz = (vp_int)(d[2]-e[2]);

    vp_int alift = ax*ax + ay*ay + az*az;
    vp_int blift = bx*bx + by*by + bz*bz;
    vp_int clift = cx*cx + cy*cy + cz*cz;
    vp_int dlift = dx*dx + dy*dy + dz*dz;

    vp_int ab = ax*by - bx*ay, bc = bx*cy - cx*by, cd = cx*dy - dx*cy;
    vp_int da = dx*ay - ax*dy, ac = ax*cy - cx*ay, bd = bx*dy - dx*by;

    vp_int abc = az*bc - bz*ac + cz*ab;
    vp_int bcd = bz*cd - cz*bd + dz*bc;
    vp_int cda = cz*da + dz*ac + az*cd;
    vp_int dab = dz*ab + az*bd + bz*da;

    {   /* Degree-5 products, accumulated at 256 bits. See vp_mul_wide. */
        vp_i256 s = vp_mul_wide(dlift, abc);
        s = vp_add256(s, vp_mul_wide(-clift, dab));
        s = vp_add256(s, vp_mul_wide(blift, cda));
        s = vp_add256(s, vp_mul_wide(-alift, bcd));
        return vp_sign256(s);
    }
}
