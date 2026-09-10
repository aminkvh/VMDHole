# Performance

Measured on an 8-core AMD Ryzen 7 7700X. Trajectory comparisons use 15
VMDPathFinder workers against serial baselines. Results depend on the system,
settings, and hardware; repetition counts and build details are recorded with
each benchmark.

| What | Compared with | Speedup |
|---|---|---|
| Surface triangulation (`sos_triangle`), same surface, byte-identical output | HOLE 2 `sos_triangle` | 3.4x at dot density 10, 73x at density 40 |
| HOLE pipeline on one structure (1BL8): search, dots, surface | HOLE 2 binaries built at -O2 | 5.3x circular, 7.1x Connolly |
| 50-frame trajectory, radius profiles | mdahole2 | 93x |
| 50-frame trajectory, radius profiles | serial HOLE shell loop | 6.1x |
| 50-frame trajectory, surface generation | mdahole2 | 50x |
| 50-frame trajectory, surface generation | serial HOLE shell loop | 16x |
| Job pool, 15 workers on 8 cores | 1 worker | 6.4x |
| Pore search step only, Nelder-Mead (18,677 atoms) | HOLE's Monte Carlo search | 5x |
| Tunnel search (MOLE 2 algorithm), matching tunnel counts | MOLE 2 | 9.7x to 16x |
| Tunnel search, compiled engine | the plugin's Tcl fallback | 173x to 231x |
| Cross-frame tunnel clustering, 1837 pathways | the plugin's Tcl fallback | 23x |

HOLE comparisons use locally rebuilt `-O2` binaries. Surface-generation timings
exclude interactive display. The Nelder-Mead 5x result measures the search
step, not a complete run; the separate 50-frame benchmark showed no overall
speedup over accelerated HOLE. Matching tunnel counts do not establish identical
route geometry; tetrahedron counts differ on some test structures.

<p align="center">
  <img src="images/performance_summary.png" alt="Benchmark summary: end-to-end trajectory throughput, triangulation, parallel scaling, tunnel search and cross-frame clustering" width="900">
</p>

[Benchmark records and reproduction instructions](https://github.com/aminkvh/VMDPathFinder/blob/main/paper/README.md).
