/* .NET Framework's seeded System.Random, ported.
 *
 * MOLE jitters every atom by +-0.00005 A from `new Random(0)` before
 * triangulating, and consumes the draws in atom order. Reproducing that
 * sequence exactly is a precondition for reproducing their triangulation, so
 * this is not an approximation of a random source - it is a bit-exact port of
 * one specific generator.
 *
 * Knuth's subtractive generator, as in the .NET reference source. Seeded
 * Random still uses it in current .NET for compatibility (only the parameterless
 * constructor moved to xoshiro), and mono implements the same, so the sequence
 * is portable.
 *
 * Verified against mono's own output for seed 0 - see mole_rng_test.c.
 */
#include "mole_rng.h"

#define MBIG  2147483647
#define MSEED  161803398

void mole_rng_init(mole_rng *r, int seed)
{
    int i, k, ii, mj, mk;
    int subtraction = (seed == (-2147483647 - 1)) ? MBIG
                    : (seed < 0 ? -seed : seed);

    mj = MSEED - subtraction;
    r->seed_array[55] = mj;
    mk = 1;

    /* Indices are walked in steps of 21 mod 55, so the array fills in a
       scattered order; index 0 is never used. */
    for (i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        r->seed_array[ii] = mk;
        mk = mj - mk;
        if (mk < 0) mk += MBIG;
        mj = r->seed_array[ii];
    }

    for (k = 1; k < 5; k++) {
        for (i = 1; i < 56; i++) {
            r->seed_array[i] -= r->seed_array[1 + (i + 30) % 55];
            if (r->seed_array[i] < 0) r->seed_array[i] += MBIG;
        }
    }

    r->inext = 0;
    r->inextp = 21;
}

static int internal_sample(mole_rng *r)
{
    int ret;
    int loc_inext = r->inext, loc_inextp = r->inextp;

    if (++loc_inext >= 56) loc_inext = 1;
    if (++loc_inextp >= 56) loc_inextp = 1;

    ret = r->seed_array[loc_inext] - r->seed_array[loc_inextp];
    if (ret == MBIG) ret--;
    if (ret < 0) ret += MBIG;

    r->seed_array[loc_inext] = ret;
    r->inext = loc_inext;
    r->inextp = loc_inextp;
    return ret;
}

double mole_rng_next_double(mole_rng *r)
{
    /* Written as the multiplication .NET uses, not a division: the two do not
       round identically and the result is compared bit-for-bit. */
    return internal_sample(r) * (1.0 / MBIG);
}

/* MOLE's general-position jitter, applied to EVERY atom in file order.
 *
 * The association is not incidental. C# evaluates
 *     a.Position.X + 0.0001 * rnd.NextDouble() - 0.00005
 * left to right, as (X + 0.0001*r) - 0.00005. Writing it as
 * `x += 0.0001*r - 0.00005` groups it as X + (0.0001*r - 0.00005) instead, and
 * floating-point addition is not associative: measured against MOLE's own
 * jittered coordinates that mis-grouping moved 3386 of 11427 values by up to
 * 7.1e-15 A. Below our grid, so it cannot move our triangulation - but it is
 * wrong, and under MOLE mode, which uses their arithmetic, it would matter. */
void mole_jitter(double *xyz, int n, int seed)
{
    mole_rng r;
    int i;
    mole_rng_init(&r, seed);
    for (i = 0; i < n; i++) {
        xyz[3*i+0] = xyz[3*i+0] + 0.0001 * mole_rng_next_double(&r) - 0.00005;
        xyz[3*i+1] = xyz[3*i+1] + 0.0001 * mole_rng_next_double(&r) - 0.00005;
        xyz[3*i+2] = xyz[3*i+2] + 0.0001 * mole_rng_next_double(&r) - 0.00005;
    }
}
