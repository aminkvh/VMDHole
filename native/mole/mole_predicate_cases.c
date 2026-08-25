/* Emit predicate signs for random and degenerate inputs, for the Tcl port to
   be checked against. Same quantised integer coordinates both sides.

   Also emits quantisation cases, prefixed "Q". Those are NOT covered by the
   predicate rows, which start from integers: vp_quant's round-half-away-from-
   zero branch is reachable only from a double whose scaled value lands exactly
   on a half step, which no real coordinate does after MOLE's jitter. Each row
   carries a flag saying whether that branch was actually taken, so the Tcl test
   can assert it was exercised rather than assume it. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../voronoi/vor_predicates.h"

static void quant_case(double x)
{
    double s = x * VP_SCALE;
    int half = (s - floor(s)) == 0.5 || (ceil(s) - s) == 0.5;
    printf("Q %.17g %ld %d\n", x, vp_quant(x), half);
}

int main(void){
    long p[15]; int i,k; unsigned s=12345;

    /* Halfway values built as (k + 0.5) / SCALE, both signs, so the tie-breaking
       branch is hit rather than hoped for; then ordinary values around it. */
    for (i = -3; i <= 3; i++) {
        quant_case((i + 0.5) / VP_SCALE);
        quant_case((i * 1000 + 0.5) / VP_SCALE);
        quant_case((double)i / VP_SCALE);
    }
    quant_case(0.0);
    quant_case(-0.0);
    for (i = 0; i < 40; i++) {
        s = s*1103515245u+12345u;
        quant_case(((double)(s % 2000001) - 1000000.0) / 8192.0);
    }

    for(i=0;i<4000;i++){
        for(k=0;k<15;k++){ s=s*1103515245u+12345u; p[k]=(long)((s>>8)%2000001)-1000000; }
        if(i%4==1){ for(k=0;k<3;k++) p[9+k]=p[k]+(p[3+k]-p[k])+(p[6+k]-p[k]); } /* coplanar */
        if(i%4==2){ for(k=0;k<3;k++) p[12+k]=p[k]; }                            /* cospherical-ish */
        printf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %d %d\n",
            p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p[9],p[10],p[11],p[12],p[13],p[14],
            vp_orient3d(p,p+3,p+6,p+9), vp_insphere(p,p+3,p+6,p+9,p+12));
    }
    return 0;
}
