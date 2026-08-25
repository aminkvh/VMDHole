/* mole_chain_cmp against MOLE's own InvariantCulture order.
 *
 * The expected string is what .NET printed under MOLE's own mono for
 * OrderBy on this list; see tests/chain_collation.tcl for the Tcl half. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mole_complex.h"

static int cmp(const void *a, const void *b)
{
    return mole_chain_cmp(*(const char *const *)a, *(const char *const *)b);
}

int main(void)
{
    const char *v[] = { "A","a","B","b","C","c","Z","z","0","1","9","?","_",
                        "AA","Ab","aB" };
    const char *want = "?,_,0,1,9,a,A,AA,aB,Ab,b,B,c,C,z,Z";
    char got[256] = "";
    int n = (int)(sizeof v / sizeof v[0]), i;
    qsort(v, (size_t)n, sizeof v[0], cmp);
    for (i = 0; i < n; i++) {
        if (i) strcat(got, ",");
        strcat(got, v[i]);
    }
    if (!strcmp(got, want)) {
        printf("  C chain collation = MOLE's InvariantCulture order      PASS\n");
        return 0;
    }
    printf("  C chain collation = MOLE's InvariantCulture order      FAIL\n");
    printf("     want %s\n     got  %s\n", want, got);
    return 1;
}
