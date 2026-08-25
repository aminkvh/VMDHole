/* .NET Framework's seeded System.Random. See mole_rng.c for why this exists. */
#ifndef MOLE_RNG_H
#define MOLE_RNG_H

typedef struct {
    int seed_array[56];   /* index 0 unused, as in the original */
    int inext, inextp;
} mole_rng;

void   mole_rng_init(mole_rng *r, int seed);
double mole_rng_next_double(mole_rng *r);

/* Apply MOLE's +-0.00005 A jitter to n atoms in place, in file order. */
void   mole_jitter(double *xyz, int n, int seed);

#endif
