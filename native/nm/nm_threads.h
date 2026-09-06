/* nm_threads.h - one thread-count policy for every tool in this directory.
 *
 * OpenMP's default is nproc, which counts hyperthreads. Measured on the pore
 * search, 18,677 atoms: 1 thread 59 ms, 4 threads 22 ms, 8 threads (this
 * machine's physical core count) 21 ms, 12 threads 24 ms, 16 threads 368 ms.
 * The last is a 17x cliff rather than a plateau - two hyperthreads on one
 * core both stream the atom grid through the same cache - so the default
 * lands on the wrong side of it. This picks physical cores instead, and
 * defers to OMP_NUM_THREADS whenever the caller sets it. Thread count never
 * changes a result: the swept outputs are byte-identical at every setting.
 */
#ifndef NM_THREADS_H
#define NM_THREADS_H
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int nm_physical_cores(void) {
    FILE *f = fopen("/proc/cpuinfo", "r");
    int siblings = 0, cores = 0, procs = 1;
    char line[256];
    if (f) {
        while (fgets(line, sizeof line, f)) {
            int v;
            if (sscanf(line, "siblings : %d", &v) == 1 && v > siblings) siblings = v;
            if (sscanf(line, "cpu cores : %d", &v) == 1 && v > cores) cores = v;
        }
        fclose(f);
    }
#ifdef _OPENMP
    procs = omp_get_num_procs();
#endif
    if (siblings > 0 && cores > 0 && siblings > cores) {
        int per = siblings / cores;
        if (per > 0 && procs / per >= 1) return procs / per;
    }
    if (cores > 0) return cores;
    return procs >= 1 ? procs : 1;
}

/* Spinning workers lose badly on a busy machine: between the many short
   per-slice parallel regions the default GOMP policy spins, and with other
   processes on the cores that measured 474 ms for a 14 ms search (HOLE under
   the same load: 1.26 s). Passive waiting gives 63 ms. The runtime reads the
   policy before main() runs, so it is set by re-executing once; if that is
   impossible the run just proceeds with the default. Call first in main(). */
#include <unistd.h>
static void nm_set_wait_policy(char **argv) {
#ifdef _OPENMP
    if (getenv("OMP_WAIT_POLICY") || getenv("NM_NO_REEXEC")) return;
    setenv("OMP_WAIT_POLICY", "PASSIVE", 1);
    setenv("NM_NO_REEXEC", "1", 1);
    execv("/proc/self/exe", argv);
    execv(argv[0], argv);
    /* neither worked: carry on with whatever the runtime chose */
#else
    (void)argv;
#endif
}

static void nm_set_default_threads(void) {
#ifdef _OPENMP
    if (getenv("OMP_NUM_THREADS")) return;      /* the caller decides */
    int n = nm_physical_cores();
    if (n >= 1) omp_set_num_threads(n);
#endif
}
#endif
