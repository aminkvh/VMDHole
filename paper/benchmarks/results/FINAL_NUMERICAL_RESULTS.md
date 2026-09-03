# FINAL NUMERICAL RESULTS - VMDHole 1.0.0 manuscript freeze

Generated: 2026-09-03T06:09:35Z
Tested commit: `07790cc70bb1d01134266e8da552e19b9a882fe7` (branch main; working tree at run time was CLEAN
except these regenerated result files themselves)
Protocol: `paper/benchmarks/reproduce.sh`, the repository's own harness,
3 timing repetitions per point, medians reported. No benchmark script,
application code, parameter, or dataset was modified for this run.

## Environment

See `env_manifest.txt` (regenerated this run) for the full record. Summary:
AMD Ryzen 7 7700X (8c/16t), 128 GB RAM, Linux 7.0.0-28-generic,
gcc/gfortran 13.3.0, Python 3.12.3, Tcl 8.6.14, VMD 2.0.1a1 (`vmd2`),
OMP_NUM_THREADS unset. Binary SHA-256s for every engine invoked are in each
CSV's own provenance header; the `env` stage verifies the stock build links
no OpenMP (the 2026-08-28 self-comparison guard).

## Commands executed (in order)

1. `bash reproduce.sh --tier1` - all 12 tier-1 stages, quiet-machine gate active
2. `bash reproduce.sh --tier2 --data vmdhole` - endtoend passed; the `scaling`
   stage's load gate then timed out at load 0.55 (>0.5 for 10 min; residual
   desktop processes) and the run stopped, per protocol
3. `bash reproduce.sh --allow-noisy --stage scaling --stage tunnel-clustering-real
   --data vmdhole` - resumed with user authorization ("noise is fine");
   loadavg at start ~0.5, recorded in the CSV headers as allow_noisy: 1
4. VMDHole-side CHAP profile regeneration per the chap stage's own staleness
   remedy (`verify_chap.tcl`), then `bash reproduce.sh --allow-noisy --stage chap`
5. `bash reproduce.sh --allow-noisy --stage figures`

Correctness gate BEFORE benchmarking: `vmdhole/tests/run_tests.sh` = ALL 20
GROUPS PASSED on `07790cc`; the tier-1 `regress` stage repeated this under
VMDHOLE_RELEASE=1 (skips forbidden): PASS (`regression.log`).

## Identity / numerical-parity gates (all PASS)

| gate | result | source |
|---|---|---|
| HOLE vs VMDHole profile identity | PASS | `profile_identity.csv` |
| Pure-Tcl fallback vs HOLE (all 3 pore methods) | PASS | `regression.log` this run (groups hole_tcl_fallback 25/25, _pore_methods, _e2e) |
| Packed-coordinate vs standard-PDB hand-off | PASS | `accel_parity.log` |
| Connolly/capsule/surface output identity | PASS | `gate_surface.log` this run + `regression.log` (capsule-incomplete) + `identity-sos` (verify.sh A/B/D/E) |
| 1BL8 pipeline: .sph md5 + triangle counts identical across stock/accel_1t/accel_nt | PASS | `pipeline_1bl8.csv` |
| Ellipse fit C-vs-Tcl parity | PASS | `regression.log` this run (ellipse-parity group) |
| Figure S2 self-agreement with CSVs | PASS (`FIG_RESULT PASS`) | `fig_performance_provenance.txt` |

STALE-FILE NOTE: `pore_tcl_fallback_parity.log`, `ellipse_parity.log`,
`hydro_qco_parity.log` and `capsule_f4.csv` in results/ are 2026-08-24
snapshot copies that no stage refreshes; the LIVE evidence for those gates
this run is `regression.log` (all 20 groups, this run's timestamp). They are
retained unmodified, per the never-replace-raw-results rule.

## End-to-end, 50-frame trajectory (Tier 2) - `endtoend.csv`

Baselines: bare HOLE = the stock `hole` binary scripted directly (serial, and
a 15-way `xargs -P15` control); mdahole2 = the MDAnalysis HOLE wrapper;
VMDHole = this plugin, 15 jobs, stock or accelerated binaries.

| deliverable | tool | median s | repetitions (s) |
|---|---|---|---|
| calc | bare_hole_serial | 5.9173 | 5.9053, 5.9173, 5.9487 |
| calc | bare_hole_parallel | 0.6894 | 0.6894, 0.6950, 0.6894 |
| calc | mdahole2 | 92.5573 | 92.1971, 95.0883, 92.5573 |
| calc | vmdhole_stock | 1.2756 | 1.2826, 1.2756, 1.2626 |
| calc | vmdhole_accel | 0.9339 | 0.9186, 0.9339, 0.9389 |
| calc | vmdhole_accel_pdb | 6.4773 | 6.4773, 6.4746, 6.7188 |
| surface | bare_hole_serial | 40.7186 | 40.7186, 40.7784, 40.4433 |
| surface | bare_hole_parallel | 4.6936 | 4.6936, 4.5999, 4.7728 |
| surface | mdahole2 | 125.2959 | 125.2959, 124.4301, 130.3889 |
| surface | vmdhole_stock | 6.2130 | 6.2347, 6.1178, 6.2130 |
| surface | vmdhole_accel | 2.2849 | 2.3236, 2.2849, 2.2769 |
| surface | vmdhole_accel_pdb | 13.5383 | 13.5355, 13.5383, 13.5395 |

**Derived ratios (full precision -> manuscript rounding):**

- calc vs mdahole2: 99.1084 -> **99.1x**
- calc vs bare serial HOLE: 6.3361 -> **6.34x**
- calc vs xargs control: 0.7382 -> **0.74x** (plugin SLOWER; manuscript's original claim stands)
- surface vs mdahole2: 54.8365 -> **54.8x**
- surface vs bare serial: 17.8207 -> **17.8x**
- surface vs xargs control: 2.0542 -> **2.05x** (plugin FASTER)
- surface accel vs plugin's own stock: 2.7192 -> **2.72x**

## Worker-count sweep (Tier 2) - `scaling.csv`

| jobs | median s | speedup vs 1 job |
|---|---|---|
| 1 | 6.4502 | 1.0000 |
| 2 | 3.4389 | 1.8757 |
| 3 | 2.4761 | 2.6050 |
| 4 | 1.9985 | 3.2275 |
| 5 | 1.6419 | 3.9285 |
| 6 | 1.5121 | 4.2657 |
| 7 | 1.3858 | 4.6545 |
| 8 | 1.2834 | 5.0259 |
| 9 | 1.1674 | 5.5253 |
| 10 | 1.0817 | 5.9630 |
| 11 | 1.0479 | 6.1554 |
| 12 | 1.0354 | 6.2297 |
| 13 | 0.9729 | 6.6299 |
| 14 | 0.9304 | 6.9327 |
| 15 | 0.9261 | 6.9649 |

15-worker headline: **0.9261 s**, **6.96x** (sub-linear; ~46% parallel efficiency).

## Triangulation-density sweep (9HNR frame 0) - `sos_scaling.csv`

| dotden | triangles | upstream ms | fast ms | speedup | identical |
|---|---|---|---|---|---|
| 10 | 2810 | 99.2 | 29.8 | 3.3 | yes |
| 15 | 6390 | 458.4 | 48.2 | 9.5 | yes |
| 20 | 11280 | 1354.2 | 79.9 | 16.9 | yes |
| 25 | 17354 | 3066.9 | 110.8 | 27.7 | yes |
| 30 | 24570 | 6055.4 | 150.6 | 40.2 | yes |
| 35 | 33388 | 11215.8 | 202.0 | 55.5 | yes |
| 40 | 43707 | 18890.2 | 262.2 | 72.0 | yes |

Range: **3.3x (density 10) to 72.0x (density 40)**, output identical at every point.

## MOLE 2 tunnel validation + timing - `tunnel_vs_mole2.csv`, `tunnel_vs_mole2_auto_origin.csv`

| structure | tetra MOLE2 | tetra VMDHole | tunnels (both) | MOLE2 s | VMDHole s | speedup |
|---|---|---|---|---|---|---|
| 1BL8 | 18394 | 18380 | 4/4 | 0.3514 | 0.0246 | 14.28x |
| 1MXT_noHET | 45044 | 45031 | 5/5 | 0.5721 | 0.0574 | 9.98x |
| 1ERI | 15163 | 15163 | 1/1 | 0.3012 | 0.0193 | 15.59x |
| 1BL8_auto | 18394 | 18380 | 8/8 | 0.3729 | 0.0282 | 13.22x |
| 1MXT_auto | 45044 | 45031 | 13/13 | 0.6716 | 0.0711 | 9.45x |

Timing range: **9.45-15.59x** (previous record 8.2-12.5x; the
tunnel engine gained an independent ~3.3x optimization merged 2026-09-01 in PR #1,
so a higher band on unchanged validation counts is the expected outcome, not drift).

## CAVER comparison - `tunnel_vs_caver_timing.csv`, `tunnel_tcl_vs_compiled.csv`, `tunnel_clustering_real_pool.csv`

| threshold | VMDHole s | CAVER stage-only s | CAVER full s | stage-only ratio | full ratio |
|---|---|---|---|---|---|
| 2.0 | 0.014005 | 0.147579 | 0.262494 | 10.54x | 18.74x |
| 4.0 | 0.010709 | 0.146122 | 0.261677 | 13.64x | 24.43x |

Stage-only band **10.5-13.6x**; complete-command band **18.7-24.4x**.

Compiled vs Tcl tunnel engine (search): 1BL8 173x; 1MXT_noHET 231x - `tunnel_tcl_vs_compiled.csv`.
Cross-frame clustering on the real 50-frame pool (1837 pathways):
compiled 0.628 s vs Tcl 16.61 s = **26.5x**,
cluster multisets identical (yes).

## CHAP replication - `chap_correlation.csv` (CHAP 0.9.1, example-02, 11 frames)

| series | metric | value |
|---|---|---|
| min_radius_per_frame | pearson_r | 0.998341 |
| min_radius_per_frame | spearman_r | 1 |
| min_radius_per_frame | rmsd | 0.0112461 |
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

Headlines (manuscript rounding): min-radius **r = 0.998**, RMSD **0.0112 A**;
registered profiles: radius **r = 0.986**, water density **r = 0.952**,
free energy **r = 0.968**.
NOTE: the run request's consistency checks listed density 0.979 / energy
0.986 - those are the pre-2026-08-24 values; the freeze record documents
their regeneration to 0.952/0.968 as the publishable numbers, and this run
reproduces exactly those.

## Per-stage pipeline (1BL8) - `pipeline_1bl8.csv`

| method | build | total s |
|---|---|---|
| circular | stock | 0.7882 |
| circular | accel_1t | 0.1491 |
| circular | accel_nt | 0.1484 |
| connolly | stock | 2.5958 |
| connolly | accel_1t | 0.7032 |
| connolly | accel_nt | 0.3605 |

circular total **5.31x**; connolly total **7.20x** (connolly at dotden 8 -
stock sos_triangle overflows above dotden ~12; always carry this caveat).

## Consistency-check reconciliation (request checklist vs this run)

| check | expected | measured | status |
|---|---|---|---|
| CHAP min radius r / RMSD | 0.998 / 0.0112 | 0.998 / 0.0112 | MATCH |
| CHAP radius r | 0.983 | 0.986 | MATCH (current record 0.986) |
| CHAP energy r | 0.986 (stale) | 0.968 | MATCHES CURRENT RECORD 0.968 |
| CHAP density r | 0.979 (stale) | 0.952 | MATCHES CURRENT RECORD 0.952 |
| 1BL8 tetrahedra / tunnels | 18394/18380, 4 | 18394/18380, 4 | EXACT |
| 1MXT tetrahedra / tunnels | 45044/45031, 5 | 45044/45031, 5 | EXACT |
| pore profiles fold | ~93.0x | 99.1x | HIGHER (see note) |
| surfaces fold | ~50.2x | 54.8x | in band |
| triangulation | 3.4-73.0x | 3.3-72.0x | in band |
| 15-worker | ~0.939 s, 6.89x | 0.9261 s, 6.96x | MATCH |
| MOLE2 fold | 8.2-12.5x | 9.4-15.6x | HIGHER - tunnel engine 3.31x fix (PR #1) |
| CAVER stage-only | 10.5-13.6x | 10.5-13.6x | EXACT |
| CAVER complete | 18.7-24.4x | 18.7-24.4x | EXACT |

Note on the pore-profile fold: mdahole2's median moved 89.4->92.6 s between
runs (its own noise) while VMDHole accel moved 0.974->0.934 s; both shifts
are within their historical run-to-run spread, and their ratio lands higher.
No code or parameter changed on either side of that comparison this cycle.

## Figure S2

`results/fig_performance.png/.pdf/.eps` regenerated by the `figures` stage
from these exact CSVs; the stage's own agreement check passed (`FIG_RESULT
PASS`, see `fig_performance_provenance.txt`, which names every source CSV
and the panel it feeds). Error bars = min/max of the 3 repetitions.
CAVER stage-only vs complete-command are separate labelled series.
