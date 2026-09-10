# FINAL NUMERICAL RESULTS - VMDPathFinder 1.0.1 freeze, Nelder-Mead search + marching-cubes mesher

Generated: 2026-09-07T01:14:34Z
Benchmarked commit: `061f7d75dafd9c4bbd3126b6ff1176844359eb49` (branch main; working tree CLEAN at run time
except these regenerated result files). The commits after it, up to `71d9ae9` where this
document was written, change only the CPOINT/CVECT stick dialog and documentation - no measured path.
On 2026-09-07 the project was renamed (commit `9ed96e9`): labels and paths in these result files
read `vmdpathfinder`/`VMDPathFinder` where the run wrote `vmdhole`/`VMDHole`. That edit is textual;
no measured value changed, and `SHA256_MANIFEST.txt` was re-hashed for it on 2026-09-10.
Protocol: `paper/benchmarks/reproduce.sh`, the repository's own harness,
3 timing repetitions per point, medians reported. No benchmark script parameter
or dataset changed for this run; the harness gained three plugin rows
(`vmdpathfinder_nm`, `vmdpathfinder_csg`, `vmdpathfinder_nm_csg`) so the plugin's own search engine
and mesher are timed beside the accelerated HOLE rows.

## Environment

See `env_manifest.txt` (regenerated this run) for the full record. Summary:
AMD Ryzen 7 7700X (8c/16t), 128 GB RAM, Linux 7.0.0-28-generic,
gcc/gfortran 13.3.0, Python 3.12.3, Tcl 8.6.14, VMD 2.0.1a1 (`vmd2`),
OMP_NUM_THREADS unset. Binary SHA-256s for every engine invoked are in each
CSV's own provenance header; the `env` stage verifies the stock build links
no OpenMP.

## Commands executed (in order)

1. `bash reproduce.sh --tier1 --tier2 --data vmdpathfinder` - env, regress, gate-surface, identity-profile, identity-sos, identity-accel and poreanalyser passed with the quiet-machine gate active; the gate then timed out at load ~2.6 (>0.5 for 10 min: the desktop's own file-sync daemon at ~2.3 cores, browser and a live VMD session, none of them benchmark processes) and the run stopped, per protocol
2. `bash reproduce.sh --tier1 --tier2 --data vmdpathfinder --allow-noisy --skip <the 7 stages above>` - the timed stages and `figures`, recorded in every CSV header as allow_noisy: 1 with the load at start
3. `bash reproduce.sh --allow-noisy --stage endtoend --stage figures --data vmdpathfinder` - the end-to-end stage again after a harness fix: VMD drops empty `-args` entries, so the optional arguments now travel as `-`. Before the fix the trailing keep-PDB flag arrived in the job-count slot, which means the `vmdpathfinder_accel_pdb` row of the 2026-09-03 freeze (6.48 s) timed a SERIAL 1-job run with the packed record, not a keep-PDB run; this run measures what the row says.
4. `bash reproduce.sh --allow-noisy --stage chap` - CHAP replication (its inputs are the checked-in CHAP 0.9.1 example-02 outputs)

Correctness gate BEFORE benchmarking: `vmdpathfinder/tests/run_tests.sh` = ALL 23
GROUPS PASSED and `tests/unit/run_unit_tests.sh` = 14 passed on `71d9ae9`; the tier-1 `regress` stage
repeated the main suite under VMDPATHFINDER_RELEASE=1: ALL 23 TEST GROUPS PASSED (`regression.log`).

## Identity / numerical-parity gates (all PASS)

| gate | result | source |
|---|---|---|
| HOLE vs VMDPathFinder profile identity | PASS (4 frames) | `profile_identity.csv` |
| Pure-Tcl fallback vs HOLE (all 3 pore methods) | PASS | `regression.log` (groups hole_tcl_fallback, _pore_methods, _e2e) |
| Packed-coordinate vs standard-PDB hand-off | PASS | `accel_parity.log` |
| Connolly/capsule/surface output identity | PASS | `gate_surface.log` + `identity-sos` (verify.sh A/B/D/E) |
| 1BL8 pipeline: .sph md5 + triangle counts identical across stock/accel_1t/accel_nt | PASS | `pipeline_1bl8.csv` |
| Ellipse fit C-vs-Tcl parity | PASS | `regression.log` (ellipse-parity group) |
| Nelder-Mead search, Connolly port, marching-cubes mesher, lobe classifier vs their references | PASS | `regression.log` (nm_engine, conn_lobes_engine, mesh_csg_engine groups) |
| Figure S2 self-agreement with CSVs | PASS (FIG_RESULT PASS) | `figures` stage output, `fig_performance_provenance.txt` |

## End-to-end, 50-frame trajectory (Tier 2) - `endtoend.csv`

Baselines: bare HOLE = the stock `hole` binary scripted directly (serial, and
a 15-way `xargs -P15` control); mdahole2 = the MDAnalysis HOLE wrapper;
VMDPathFinder = this plugin, 15 jobs, stock or accelerated binaries. The three new
rows keep the accelerated binaries and switch the plugin's own engines on:
`vmdpathfinder_nm` = Nelder-Mead search instead of HOLE's Monte Carlo (calc);
`vmdpathfinder_csg` = marching-cubes mesher on HOLE's search (surface);
`vmdpathfinder_nm_csg` = both (surface).

| deliverable | tool | median s | repetitions (s) |
|---|---|---|---|
| calc | bare_hole_serial | 6.1398 | 6.1611, 6.1346, 6.1398 |
| calc | bare_hole_parallel | 0.7987 | 0.7987, 0.7973, 0.8117 |
| calc | mdahole2 | 93.8384 | 94.2621, 93.8384, 93.8282 |
| calc | vmdpathfinder_stock | 1.4582 | 1.4562, 1.4582, 1.4723 |
| calc | vmdpathfinder_accel | 1.0047 | 0.9715, 1.0188, 1.0047 |
| calc | vmdpathfinder_accel_pdb | 1.4255 | 1.4672, 1.4175, 1.4255 |
| calc | vmdpathfinder_nm | 1.0269 | 1.0144, 1.0269, 1.0298 |
| surface | bare_hole_serial | 42.6641 | 42.6641, 41.9096, 42.7141 |
| surface | bare_hole_parallel | 5.2479 | 5.1730, 5.2575, 5.2479 |
| surface | mdahole2 | 130.9052 | 130.9052, 130.7011, 131.6016 |
| surface | vmdpathfinder_stock | 6.7874 | 6.7874, 6.7734, 6.8478 |
| surface | vmdpathfinder_accel | 2.6224 | 2.6676, 2.6224, 2.5101 |
| surface | vmdpathfinder_accel_pdb | 2.8904 | 2.9487, 2.8904, 2.8691 |
| surface | vmdpathfinder_csg | 2.0426 | 2.1177, 2.0426, 2.0033 |
| surface | vmdpathfinder_nm_csg | 2.0750 | 2.1006, 2.0630, 2.0750 |

**Derived ratios (full precision -> rounding):**

- calc vs mdahole2 (accelerated HOLE search): 93.3994 -> **93.4x**
- calc vs mdahole2 (Nelder-Mead search): 91.3803 -> **91.4x**
- calc, Nelder-Mead vs accelerated HOLE search: 0.9784 -> **0.98x**
- calc vs bare serial HOLE: 6.1111 -> **6.11x**
- calc vs xargs control: 0.7950 -> **0.79x**
- surface vs mdahole2 (accelerated sph_process + sos_triangle): 49.9181 -> **49.9x**
- surface vs mdahole2 (Nelder-Mead + marching cubes): 63.0868 -> **63.1x**
- surface, marching cubes vs sos_triangle on the same HOLE search: 1.2839 -> **1.28x**
- surface, Nelder-Mead + marching cubes vs accelerated HOLE path: 1.2638 -> **1.26x**
- surface vs bare serial: 16.2691 -> **16.27x**
- surface vs xargs control: 2.0012 -> **2.00x**
- surface accel vs plugin's own stock: 2.5882 -> **2.59x**

## Worker-count sweep (Tier 2) - `scaling.csv`

| jobs | median s | speedup vs 1 job |
|---|---|---|
| 1 | 6.4691 | 1.0000 |
| 2 | 3.4569 | 1.8714 |
| 3 | 2.4861 | 2.6021 |
| 4 | 2.0010 | 3.2329 |
| 5 | 1.6799 | 3.8509 |
| 6 | 1.6736 | 3.8654 |
| 7 | 1.5773 | 4.1014 |
| 8 | 1.3395 | 4.8295 |
| 9 | 1.2086 | 5.3526 |
| 10 | 1.1281 | 5.7345 |
| 11 | 1.0990 | 5.8864 |
| 12 | 1.2111 | 5.3415 |
| 13 | 1.0482 | 6.1716 |
| 14 | 1.0465 | 6.1817 |
| 15 | 1.0123 | 6.3905 |

15-worker headline: **1.0123 s**, **6.39x**.

## Triangulation-density sweep (9HNR frame 0) - `sos_scaling.csv`

| dotden | triangles | upstream ms | fast ms | speedup | identical |
|---|---|---|---|---|---|
| 10 | 2810 | 105.5 | 30.9 | 3.4 | yes |
| 15 | 6390 | 464.2 | 50.5 | 9.2 | yes |
| 20 | 11280 | 1399.1 | 79.0 | 17.7 | yes |
| 25 | 17354 | 3216.1 | 114.3 | 28.1 | yes |
| 30 | 24570 | 6400.7 | 159.9 | 40.0 | yes |
| 35 | 33388 | 11644.1 | 208.3 | 55.9 | yes |
| 40 | 43707 | 19775.4 | 270.1 | 73.2 | yes |

Range: **3.4x (density 10) to 73.2x (density 40)**, output identical at every point.

## MOLE 2 tunnel validation + timing - `tunnel_vs_mole2.csv`, `tunnel_vs_mole2_auto_origin.csv`

| structure | tetra MOLE2 | tetra VMDPathFinder | tunnels (both) | MOLE2 s | VMDPathFinder s | speedup |
|---|---|---|---|---|---|---|
| 1BL8 | 18394 | 18380 | 4/4 | 0.3813 | 0.0250 | 15.24x |
| 1MXT_noHET | 45044 | 45031 | 5/5 | 0.6124 | 0.0581 | 10.55x |
| 1ERI | 15163 | 15163 | 1/1 | 0.3224 | 0.0197 | 16.41x |
| 1BL8_auto | 18394 | 18380 | 8/8 | 0.3979 | 0.0287 | 13.87x |
| 1MXT_auto | 45044 | 45031 | 13/13 | 0.7249 | 0.0750 | 9.66x |

Timing range: **9.66-16.41x**; tunnel counts agree on every structure: yes.

## CAVER comparison - `tunnel_vs_caver_timing.csv`, `tunnel_tcl_vs_compiled.csv`, `tunnel_clustering_real_pool.csv`

| threshold | VMDPathFinder s | CAVER stage-only s | CAVER full s | stage-only ratio | full ratio |
|---|---|---|---|---|---|
| 2.0 | 0.014005 | 0.147579 | 0.262494 | 10.54x | 18.74x |
| 4.0 | 0.010709 | 0.146122 | 0.261677 | 13.64x | 24.43x |

Stage-only band **10.5-13.6x**; complete-command band **18.7-24.4x**.

Compiled vs Tcl tunnel engine (search): 1BL8 173.4x; 1MXT_noHET 231.4x - `tunnel_tcl_vs_compiled.csv`.
Cross-frame clustering on the real 50-frame pool (1837 pathways), the kernel now answered by the resident mesher process:
compiled 0.774 s vs Tcl 18.02 s = **23.3x**,
cluster multisets identical (yes).

## CHAP replication - `chap_correlation.csv` (CHAP 0.9.1, example-02, 11 frames)

| series | metric | value |
|---|---|---|
| min_radius_per_frame | pearson_r | 0.998341 |
| min_radius_per_frame | spearman_r | 1 |
| min_radius_per_frame | rmsd | 0.0112461 |
| min_radius_per_frame | mean_ours | 2.06829 |
| min_radius_per_frame | mean_chap | 2.06834 |
| min_radius_per_frame | frames | none |
| registration | sigma | -1 |
| registration | delta | -3 |
| registration | radius_pearson_r | 0.985819 |
| radius | pearson_r | 0.985819 |
| radius | spearman_r | 0.980819 |
| radius | rmsd | 0.140375 |
| density | pearson_r | 0.952111 |
| density | spearman_r | 0.944332 |
| density | rmsd | 0.146735 |
| energy | pearson_r | 0.96833 |
| energy | spearman_r | 0.953239 |
| energy | rmsd | 0.226769 |
| density2 | pearson_r | 0.923326 |
| density2 | spearman_r | 0.889271 |
| density2 | rmsd | 0.197276 |
| energy2 | pearson_r | 0.955671 |
| energy2 | spearman_r | 0.902834 |
| energy2 | rmsd | 1.35709 |

Headlines: min-radius **r = 0.998**, RMSD **0.0112 A**;
registered profiles: radius **r = 0.986**, water density **r = 0.952**,
free energy **r = 0.968**.

## Per-stage pipeline (1BL8) - `pipeline_1bl8.csv`

| method | build | total s |
|---|---|---|
| circular | stock | 0.8273 |
| circular | accel_1t | 0.1575 |
| circular | accel_nt | 0.1566 |
| connolly | stock | 2.7219 |
| connolly | accel_1t | 0.7212 |
| connolly | accel_nt | 0.3857 |

circular total **5.28x**; connolly total **7.06x** (connolly at dotden 8 -
stock sos_triangle overflows above dotden ~12; always carry this caveat).

## Comparison with the 2026-09-03 freeze (commit 07790cc)

| headline | 2026-09-03 | this run |
|---|---|---|
| calc vs mdahole2 (accelerated HOLE search) | 99.1 | 93.40 |
| surfaces vs mdahole2 (accelerated HOLE path) | 54.8 | 49.92 |
| 15-worker speedup | 6.96 | 6.39 |
| triangulation, density 40 | 72.0 | 73.20 |
| MOLE 2 fold, 1BL8 | 14.28 | 15.24 |
| CAVER stage-only, 2.0 A | 10.54 | 10.54 |
| cross-frame clustering compiled vs Tcl | 26.5 | 23.27 |
| CHAP min-radius r | 0.998 | 1.00 |

The reference tools (mdahole2, bare HOLE, MOLE 2, CAVER, CHAP) and the stock
build are unchanged; movements in those rows are run-to-run spread. New this
run: the `vmdpathfinder_nm`, `vmdpathfinder_csg` and `vmdpathfinder_nm_csg` rows above and their
bars in Figure S2.

## Figure S2

`results/fig_performance.png/.pdf/.eps` regenerated by the `figures` stage
from these exact CSVs; the stage's own agreement check passed (`FIG_RESULT
PASS`, see `fig_performance_provenance.txt`). `docs/images/
performance_summary.png` is a copy of the PNG.

Nelder-Mead search and marching-cubes mesher, 15-worker pool (`endtoend.csv`):
search 0.98x the accelerated HOLE search; mesher 1.28x the accelerated
`sph_process` + `sos_triangle` path; both 1.26x. Not in Figure S2.
