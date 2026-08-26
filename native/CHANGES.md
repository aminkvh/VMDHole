# Changes relative to upstream HOLE2

This is the Apache-2.0 statement of modifications for everything in
`native/`. It has two parts: the **C program** `sos_triangle_fast.c`
(sections below, up to "Additive feature: batched recolour"), and the **Fortran
patches** in `connolly_patches/`, which modify `hole` and `sph_process` (section
"Fortran patches", at the end). Every modification in both parts is intended to
preserve output exactly; where that guarantee is qualified, the qualification is
stated with the change.

## The C program: `sos_triangle_fast.c`

`sos_triangle_fast.c` is a modified copy of `sos_triangle.c` from HOLE2
(https://github.com/osmart/hole2, Apache-2.0). It is a **drop-in replacement**:
the command-line interface, the stdin `.sos` input format, and the stdout
`.vmd_plot` output are unchanged, and the generated surface is **byte-for-byte
identical** to the upstream program (see `benchmarks/` and `verify.sh`). The
modifications are confined to performance and capacity; the triangulation
algorithm and its output are not altered.

This document describes each change for a reader/reviewer. Line-level rationale
is also in the source comments.

> **Additive feature (opt-in):** the binary also gains a hydrophobicity-colouring
> mode used by the VMDHole plugin (`--hydro-atoms`, see the last section). It is
> **inert unless that flag is passed** — with the normal `-s`/`-v` invocation the
> program is still byte-for-byte the upstream output, so the drop-in guarantee
> above is unaffected.

## Background: where the time went

`sos_triangle` turns a cloud of surface dots (produced by `sph_process`) into a
triangulated mesh using an incremental, advancing-front triangulation. Profiling
the unmodified program on a typical surface (HOLE dot density 15, ~5 000 dots)
showed the runtime concentrated in two routines, both effectively **O(N²)** in
the number of surface dots `N`:

| routine | role | share of runtime | why it is O(N²) |
|---|---|---|---|
| `neighbour()` | for each advancing edge, find the best next vertex | ~83 % | scans **all** dots for **every** edge |
| `destroy()` | mark an edge as already-triangulated | ~14 % | scans the **whole** append-only edge list for **every** lookup |

The upstream source itself flags the first one: *"search through all dots for
neighbours. (slow, but works...) Box type data structure would be better..."*.
The two superficially-quadratic routines `cull_coords()` and `check_point()`
turned out to be negligible at `-O2`, but were tidied anyway.

## The six changes

| # | Routine | Before | After | Complexity | Practical effect |
|---|---------|--------|-------|------------|------------------|
| 1 | `cull_coords` | exhaustive pairwise duplicate test | spatial hash grid (cell = the 1e-3 dedup tolerance; exact 3×3×3 neighbour test) | O(N²) → ~O(N) | negligible |
| 2 | `check_point` | rescans every triangle to ask "is this dot used yet?" | O(1) `point_used[]` flag set when a vertex is written in `gen_triangle` | O(N·T) → O(1) | negligible |
| 3 | `neighbour` (micro) | recomputes edge-only terms per candidate; uses `pow(v,2)` | hoist edge-invariant terms out of the loop; `pow(v,2)` → `v*v` | constant-factor | ~2.2× |
| 4 | `destroy` | linear scan of the whole edge list | hash keyed on the (unordered) endpoint pair; bucket chains kept in creation order so the first match is unchanged | O(E) → ~O(1) | the largest single win at typical sizes |
| 5 | `neighbour` (grid) | scans all N dots | queries a uniform spatial grid for the only dots that can win (within 3·\|base\| of vertex `a`) | O(N) per call → ~O(1) | the dominant win on large surfaces |
| 5b | `neighbour` (grid **cell size**) | cell = `cbrt(bbox_volume/N)`, which assumes 3-D filling and overestimates the 2-D surface spacing ~5× | cell = ~5× the median nearest-neighbour distance (the true surface spacing) | constant-factor, but a large one | **~2× on the triangulation phase**; the dominant win at high density (see below) |
| 6 | `MAX_COORD` | `30000` (fixed array bound) | `200000` | capacity, not speed | large pores/high densities that previously aborted now complete |

Changes 3–5b are the speedups; 1–2 are clean-ups; 6 is a capacity increase that
is independent of the speed work (it removes the "Maximum number of polygons
exceeded" abort).

> **Why 5b matters (profiled 2026-06).** After changes 1–5, `gprof`/`perf` and
> per-phase timing put the entire remaining cost in the advancing-front loop,
> and inside it `nb_consider()` was called ~220 M times on a 42 k-triangle
> surface. Instrumenting the grid query showed **>99 % of those candidates were
> outside the 3·\|base\| radius** — pure waste — because the grid cells were ~5×
> too coarse: `cbrt(volume/N)` treats the dots as filling a 3-D box, but they lie
> on a 2-D shell. Sizing cells to the true surface spacing (≈5× the median
> nearest-neighbour distance, validated near-optimal across surfaces and
> proteins) cuts the candidates ~3× and the triangulation phase ~2×. It is a pure
> acceleration: `nb_consider()` applies the identical 3·\|base\| test regardless
> of cell size, and a mis-estimate falls back to the original full scan, so the
> output is byte-for-byte unchanged (re-verified on all fixtures).

### Why the output is provably identical

Two of the changes reorder work, so they need an argument that the result cannot
change:

* **Change 5 (spatial grid) is lossless by construction.** The original
  `neighbour()` already rejects any candidate whose distance from edge-vertex
  `a` exceeds `3·|base|` (the `check_dist` test). A dot outside that radius can
  therefore never be selected and never affects the running minimum, so
  restricting the search to the grid cells covering that radius removes only
  work that had no effect. Because dots are now visited in cell order rather
  than ascending-index order, the selection picks the smallest subtended angle
  and **breaks ties by smallest dot index** — exactly the choice the original
  strict-`<` scan in index order made. A guard reproduces the original behaviour
  that the initial angle bound (1.0) is only ever beaten by a strictly smaller
  angle, never matched.

* **Change 4 (edge hash) preserves first-match order.** Edges are appended to
  their hash bucket in creation order, and the lookup returns the first matching
  endpoint pair in that bucket — the same node the original full-list scan
  returned, so the flags it sets on the edge tree are identical.

Changes 1–3 and 6 are value-preserving by inspection (same arithmetic, fewer
times; `pow(x,2)` and `x*x` agree to the printed precision and were verified
identical).

## Two speedup figures: component vs. end-to-end

It is important not to conflate these.

* **Component speedup (`sos_triangle` alone).** The triangulation program itself
  runs **~9.4× faster at the plugin's default dot density (15) and up to ~75× at
  density 40** (vs the unmodified upstream source), the speedup growing with
  surface size. Measured on **9HNR**, all rows byte-identical. These figures are
  transcribed from `../paper/benchmarks/results/sos_scaling.csv`, which is the
  data of record — if the two ever disagree, the CSV wins:

  | dot density | triangles | upstream | fast | speedup |
  |---:|---:|---:|---:|---:|
  | 10 | 2 869 | 99.9 ms | 31.6 ms | 3.2× |
  | 15 | 6 514 | 444.2 ms | 47.3 ms | 9.4× |
  | 20 | 11 452 | 1 333.5 ms | 75.1 ms | 17.8× |
  | 25 | 17 270 | 3 030.2 ms | 106.1 ms | 28.6× |
  | 30 | 24 841 | 6 180.7 ms | 147.3 ms | 42.0× |
  | 35 | 33 612 | 11 214.5 ms | 199.2 ms | 56.3× |
  | 40 | 43 712 | 18 687.4 ms | 248.6 ms | **75.2×** |

  > **Timings are medians of three and move by a few percent between runs**; the
  > table above is re-transcribed whenever the sweep is re-run, and the CSV is
  > always the tie-breaker. The **triangle counts do not move** — they are
  > byte-reproducible, which is what makes the ratios comparable at all. An
  > earlier run of this same sweep gave 9.5× and 74.5× at densities 15 and 40
  > against 9.4× and 75.2× here; treat the leading digit, not the decimal, as
  > the claim.

  > These numbers replace an earlier sweep measured on a GABA-A `.sph` taken from
  > a trajectory too large to distribute. The input here is **9HNR frame 0**, a
  > public structure, and the table is regenerated by
  > `paper/benchmarks/reproduce.sh --stage sos-sweep` rather than transcribed by hand. Triangle
  > counts are therefore lower, and the default-density figure moved
  > 13.6× → 9.4× while density 40 moved 71.3× → 75.2×. Same code, different input.
  >
  > ⚠️ **The benchmark inputs are not in the repository yet**: `9HNR.pdb` is
  > untracked and `hole_output_9HNR/` is `.gitignore`d as generated output, so a
  > clean clone cannot run this stage today. See "Benchmark inputs are not
  > committed" in `vmdhole/HANDOFF.md`. That is a packaging gap to close before
  > submission, not a property of the measurement.

* **End-to-end speedup (per analyzed frame).** Inside VMDHole each frame runs a
  three-stage pipeline — `hole` (channel search) → `sph_process` (surface
  sampling) → `sos_triangle` (triangulation).

  > ⚠️ **The table below is historical.** It was measured when `sos_triangle` was
  > the *only* accelerated stage, and is kept because it is what motivates the
  > Amdahl argument that follows. It does **not** describe the shipped binaries:
  > `hole` and `sph_process` have since been accelerated too. For the measured
  > whole-pipeline figures, see "End-to-end, re-measured" below.

  On a representative frame (~5 000 surface dots, dot density 15), with only
  `sos_triangle` replaced:

  | stage | before (stock sos) | after (fast sos) |
  |---|---:|---:|
  | `hole` | 229 ms (25%) | 229 ms (70%) |
  | `sph_process` | 27 ms (3%) | 27 ms (9%) |
  | `sos_triangle` | 675 ms (72%) | 70 ms (21%) |
  | **total per frame** | **931 ms** | **326 ms** |

  So the **end-to-end per-frame speedup is ~2.9×** (931 → 326 ms). By Amdahl's
  law this is the expected result of a ~10× speedup applied to the ~72% of the
  work that `sos_triangle` represented: `1 / (0.28 + 0.72/10) ≈ 2.8×`. After the
  change the per-frame cost is dominated by `hole`'s channel search (now ~70%),
  which is unmodified; `sos_triangle` is no longer the bottleneck. Pushing
  `sos_triangle` further therefore has limited end-to-end benefit at typical
  sizes (though it continues to help large/high-density surfaces, where the
  triangulation share is larger).

### End-to-end, re-measured against genuine stock

Stock = unpatched HOLE 2.2 built at the **same `-O2`** (`stock_build/`), 1 thread.
Accelerated = the shipped binaries, `OMP_NUM_THREADS=8`. 1BL8, seeded (`raseed 1`),
same machine, output verified identical. Raw data:
`../paper/benchmarks/pipeline_end_to_end.csv`.

| mode | dotden | stage | stock | accel | speedup |
|---|---:|---|---:|---:|---:|
| circular | 15 | `hole` | 38 ms | 37 ms | 1.0x |
| | | `sph_process` | 17 ms | 18 ms | 1.0x |
| | | `sos_triangle` | 2 402 ms | 153 ms | 15.7x |
| | | **total** | **2 457 ms** | **208 ms** | **11.8x** |
| connolly | 8 | `hole` | 103 ms | 82 ms | 1.3x |
| | | `sph_process` | 69 ms | 33 ms | 2.1x |
| | | `sos_triangle` | 1 186 ms | 82 ms | 14.5x |
| | | **total** | **1 358 ms** | **197 ms** | **6.9x** |

Two things this makes explicit that the old ~2.9x figure did not:

* **The end-to-end gain depends on what dominates.** On a small circular pore the
  surface step is ~98% of the work, so the total tracks `sos_triangle` closely
  (11.8x). Under CONNOLLY the channel search and dot generation matter more, so
  the total is lower (6.9x) even though every stage got faster — the Fortran
  parallelism is what keeps `hole` and `sph_process` from becoming the new
  bottleneck.
* ⚠️ **CONNOLLY is measured at dotden 8, not the plugin default of 15, because
  stock cannot complete it at 15.** Upstream `sos_triangle` has `MAX_COORD 30000`
  and emits **zero triangles** from dotden 12 upward on this structure (measured:
  3 714 triangles at dotden 4, 13 417 at 8, **0 at 12**) — it hits the polygon cap
  and produces no surface at all. Above that density there is no stock result to
  divide by; the difference is capability, not speed. This is why the plugin
  derives its dot budget from the binary's capabilities rather than assuming the
  fork's limits.

## Scalability and parallelism

* **Per-surface scalability.** The change is algorithmic, not just a constant
  factor: a log–log fit of runtime vs. triangle count over all seven rows of
  `../paper/benchmarks/results/sos_scaling.csv` gives exponents of **1.93
  (upstream) → 0.78 (this work)** (reproduce with
  `paper/benchmarks/reproduce.sh --stage sos-sweep`).

  A plain log–log fit understates the true exponent here, because it charges
  fixed per-run cost (process start, file I/O) to the smallest surfaces. Fitting
  the model that separates it, `t = c + a·N^b`, over the same seven rows
  resolves both terms (R² > 0.9998):

  | build | fixed cost `c` | exponent `b` |
  |---|---:|---:|
  | upstream | not resolved (see below) | **1.95** (near-quadratic) |
  | this work | 15.0 ms | **1.02** (near-linear) |

  The **exponents are the stable quantity** — across re-runs of this sweep they
  move only in the third digit (upstream 1.95–1.97, this work 1.02) while the
  timings themselves move a few percent. The fixed-cost term is well determined
  for the fast build (~13–15 ms) but **not** for upstream, whose runtime is so
  dominated by the quadratic term that the fit puts `c` near zero and it can come
  out slightly negative; quote the exponent, not upstream's `c`.

  So the change really is algorithmic — near-quadratic → near-linear — and the
  0.78 above is an artefact of the offset, not sub-linear work. Consequently the
  component speedup *grows* with surface size — 9.4× at the plugin's default
  density (15) and 75.2× on a ~44 000-triangle surface — instead of being a
  fixed multiplier.

* **Across-frame parallelism (already in the plugin).** A single `sos_triangle`
  invocation is inherently sequential (the advancing front depends on the dots
  already placed), so it is not internally multi-threaded. However, a trajectory
  produces one independent surface *per frame*, and the VMDHole plugin already
  builds these in parallel across CPU cores (`prebuild_surfaces_parallel` /
  `run_shell_pool`). This is the embarrassingly-parallel part of the workload
  and scales near-linearly with cores; the optimisation here reduces the
  *per-frame* serial cost that those workers each pay. The two are complementary:
  faster single-surface build × N-core frame parallelism.

## Verification

Output was confirmed byte-for-byte identical to the unmodified upstream program
on: the 50 real `.sos` surfaces under `vmdhole/hole_output/`; HOLE's own
example structures (gramicidin 1grm, cholera toxin 1chb) run through the full
`hole → sph_process → sos_triangle` pipeline; and dot densities 5–30. For
surfaces where the as-shipped upstream binary overflows its 30000-polygon cap
(e.g. cholera toxin at density 15), the comparison was made against the upstream
**source** rebuilt with only `MAX_COORD` raised — still identical. Re-run it
yourself with `./verify.sh`.

## Additive feature: hydrophobicity colouring (opt-in)

The VMDHole plugin can colour the pore surface by the hydrophobicity of the
nearest channel-lining residues (Kyte–Doolittle or Wimley–White). It originally
did this in a pure-Tcl post-pass (`colorize_hydrophobic`) that re-read the mesh
and, for every triangle, searched the channel spheres — fine for one frame but a
visible per-frame cost when scrubbing or playing a trajectory, and it runs on
VMD's single main thread.

This binary can now do that colouring **inside the same compiled pass that
triangulates the surface**, so it also rides the plugin's existing across-frame
parallelism (`prebuild_surfaces_parallel`). The plugin auto-detects the feature
and falls back to its Tcl path when the binary doesn't have it, so
hydrophobicity never *depends* on this binary — only its speed does.

### Interface

| flag | meaning |
|---|---|
| `--hole-features` | Print `hole_features: hydro` to stdout and exit. The plugin's capability probe; a stock binary prints help instead, so it reads as "not capable". |
| `--hydro-atoms FILE` | Enable colouring. `FILE` has one row per channel-local atom: `x y z h_kd h_ww` (positions plus the atom's hydropathy on each scale). The plugin dumps this with an `atomselect` so the scales stay single-sourced in Tcl. |
| `--hydro-sph FILE` | The HOLE `.sph` file (channel sphere centres + radii). |
| `--hydro-scheme kd\|ww` | Which hydropathy column to use (default `kd`). |
| `--points` | Emit the surface as unique vertices (`draw point`) instead of triangles — the plugin's `dots` display. |

With no `--hydro-atoms` / `--points`, none of this code runs and the output is
unchanged.

`--points` replaces the plugin's Tcl `dots_from_trinorm` pass (which read the
triangulated mesh back and de-duplicated its vertices, ~50 ms/frame on VMD's main
thread) with a direct emit from the already-de-duplicated dot list — so the dot
surface is produced in the same compiled, parallel build as the mesh. Output is
identical to `dots_from_trinorm` (same points, colours and order; `verify.sh`
Part D checks `--points` emits exactly the trinorm mesh's unique vertices).

### What it computes (identical to the Tcl path)

For each sphere, the mean hydropathy of atoms within `r + 3 Å`; the spheres are
then thinned to ≤200 (keeping the last). Each emitted triangle is coloured by
its centroid's nearest thinned sphere, mapped to a VMD colour name by the same
thresholds as `::VMDHole::hydro_to_vmd_color`. The triangle centroid is taken
from the coordinates **rounded to the 3 decimals the mesh is written with**, so
the result is *bit-identical* to the Tcl path (which reads vertices back from the
printed mesh), not merely close.

### Verification

`verify.sh` Part C runs the full `sph_process → sos_triangle_fast --hydro`
pipeline against `hydro_reference.py` (a faithful port of the Tcl colouring) on
the committed `.sph` surfaces. Measured result: **0 differing triangles out of
~85 000** across four frames × both scales.

## Additive feature: batched recolour (`--batch-recolor`)

The plugin's "Pre-build all surfaces after a run" option already batches the
triangulation phase (`--batch`, one process per worker instead of one per
frame, round-robin distributed). Its **property/hydrophobicity recolour**
phase did not have the same treatment: it recoloured an already-triangulated
base mesh (`--recolor` + `--hydro-values`) with one process **per frame**,
even though each individual recolour is fast — no triangulation, just a
nearest-sphere colour lookup per triangle centroid. At that per-call cost,
process-spawn overhead (fork/exec/shell/dynamic-link) is a proportionally
larger fraction of the total time than it is for the (much heavier)
triangulation phase, so batching it matters even though the per-job work is
cheap.

`--batch-recolor FILE` runs multiple recolour jobs in one process, values-path
only (the CLI's legacy `--hydro-atoms` atom-averaging path is not supported in
batch mode — the plugin only reaches for this feature via the values path
already). Each line: `base_vmd_plot<TAB>out_vmd_plot<TAB>values_path<TAB>sph_path`.
`--hydro-signed` applies to every job in the file, parsed before `--batch-recolor`
consumes the rest of argv (same convention as `--batch`'s global flags).

### Verification

Tested against real generated surfaces (`vmdhole/hole_output_step5_assembly.hmr/`,
three frames, ~1450 spheres / ~17-20k mesh lines each): ran each frame through
the existing single-job `--recolor --hydro-values` path, then through
`--batch-recolor` with the same three jobs **listed out of generation order**
(specifically to catch any state leaking from one job into the next via
`reset_hydro_state()` not being called, or being incomplete). All three
outputs were byte-for-byte identical to the single-job reference.

## Robustness fixes (no effect on valid input)

Three defensive changes. None alters behaviour on well-formed input; each was
verified byte-for-byte identical against the previous binary on real `.sos` data
before deployment.

### 1. `calc_tri` recursion: dedicated stack + depth ceiling

`calc_tri` recurses once per triangle, so on a large surface the recursion depth
tracks the mesh size. On the default 8 MB thread stack a big Connolly cloud
could overflow it, which manifests as a crash with no diagnostic. The recursion
is now entered through `calc_tri_root()`, which runs it on a dedicated pthread
with a 768 MB stack (`CT_STACK_BYTES`), falling back to a direct call if the
thread cannot be created. A `CT_MAX_DEPTH` counter (4,000,000 — about 30x the
deepest legitimate run measured here) turns a runaway recursion into a reported
error instead of a stack smash. The triangulation itself is untouched.

### 2. `-X` argument arity

The short-option `switch` reads the operand of `-X` from `argv[2]`. Unlike the
`--` long options — which are guarded by the `one_arg`/`two_arg` tables — it had
no arity check, so a bare `-X` with no following argument dereferenced past the
end of `argv` and segfaulted (confirmed under AddressSanitizer). An `argc < 3`
check now reports the missing operand and exits 1. `argv[1]` is the current
option and `argv[2]` its operand, so `argc >= 3` is the correct condition;
verified that `-X 5.0` still parses in any position and in combination with
other options. `-X` is the only short option that takes an operand.

### 3. `read_cord` honours the `sscanf` field count

`read_cord` ignored the return value of its `sscanf`. Because `sscanf` leaves
unmatched variables **unchanged**, a short or malformed record inherited
uninitialised stack values on the first line and the *previous* line's values on
every line after — producing geometry from data that was never in the file. The
return value is now captured and each record type requires the fields it
actually consumes: a colour record (type 1) needs 2, a point (type 4) needs 4
(type + x,y,z).

Note the threshold is 4, not 7: requiring all seven would discard records whose
geometry is perfectly usable.

> **Correction.** An earlier version of this note said fields 5-7 are
> "unconditionally overwritten with normals recomputed from the finished
> triangulation, so they are vestigial". They are not. The recomputation is
> unreachable - `vertex_normals()` returns before the loop that assigns
> `dots[][4]` (that dead code is verbatim upstream) - so the normal read from
> the file flows through to the output, and `reorder_triangle()` sums the three
> vertex normals and swaps two corners when the dot product is negative. A short
> point record inheriting the previous record's normal could therefore change
> emitted vertex order. `read_cord` now zeroes fields 5-7 before each scan, so a
> short record reads as "no normal" rather than "the last one's". A well-formed
> `.sos` always supplies seven fields, so valid input is unaffected.

## Additive feature: ellipse-probe asymmetry performance (`--asym-ellipse`)

Reported symptom: "the C accelerated ellipse is still slow — for a 10K
[trajectory] analysis it would take 10 hours" (~3.6 s/frame).

### Measurement first

Profiled a real 508-slice, 27,315-atom channel (`traj50_bare_inputs/frame_00000`,
the paper benchmark fixture) end to end with manual phase timers
(`SOS_ELLIPSE_TIMING=1` env var, stderr only, see `now_sec()`/`g_time_gather`/
`g_time_fit`) plus `/usr/bin/time`, `--asym-threads 1` (no parallelism, so wall
time = phase cost, not thread-seconds):

| phase | time | share |
|---|---|---|
| parse (`hydro_read_spheres` + `hydro3d_read_atoms_lining`) | 0.008 s | 0.2% |
| neighbour gather (per-slice atom-window search) | 0.033 s | 0.9% |
| PoreAnalyser Nelder-Mead fit (`pa_fit`/`sci_nm`/`sci_pen`) | 3.85 s | 98.9% |

So the earlier "distance matrix vs Tcl marshalling" and "process-per-frame
overhead" failure modes documented elsewhere in this file are **not** the cause
here: parsing and the neighbour search are both negligible. Instrumentation
added at this step (call counters on `sci_pen`/`pa_ev`) confirmed the fit issues
~26M `pa_ev` calls per frame (508 slices × ~57 atoms/window × ~900 penalty
evaluations/slice from 3 Nelder-Mead runs × ~127 iterations/run), of which
essentially all take the "off-ellipse" branch into `pa_de` (PoreAnalyser's
3-iteration Chatfield closest-point solver — genuine required work, left
untouched: it is the algorithm, not overhead). The counters were removed again
before shipping (see "value-preserving" below); the numbers above are what they
measured.

### Two changes, both provably value-preserving

1. **Hoist the ellipse rotation and rotated centre out of the per-atom loop.**
   `pa_on`/`pa_ev` used to recompute `cos(-th)`, `sin(-th)`, `cos(th)`, `sin(th)`
   and the rotated centre `(cx·c_nth-cy·s_nth, cx·s_nth+cy·c_nth)` on **every**
   atom, even though theta and the centre are fixed for the whole atom loop in
   both call sites (`sci_pen`: one Nelder-Mead penalty evaluation; `pa_fit`'s
   probe-shrink loop: theta is the literal `0.0`, centre is `(px,py)`, for the
   *entire* shrink loop, not just one atom pass). Each caller now computes
   `c_nth,s_nth,c_th,s_th,rc1,rc2` **once** via the exact same expressions, and
   `pa_on`/`pa_ev` take them as parameters instead of recomputing internally.
   `cos()`/`sin()`/multiply/add are pure deterministic functions of a
   bit-identical input, so reusing one evaluation is bit-for-bit identical to
   repeating it — not an approximation, a de-duplication.
2. **Auto-detect thread count for the single-frame CLI paths only.**
   `--asym-threads` already existed (OpenMP over the per-slice loop) but had to
   be passed explicitly by the caller; if the flag is never given at all (a
   stale capability-probe cache reporting `asymthreads` unsupported, or an
   older Tcl side that predates the flag), `--asym-ellipse`/`--asym-ellipse-geo`
   ran fully serial with no way to opt back in. `asym_threads_explicit`
   distinguishes "the flag was given" from "still at the built-in default":
   when unset, the two **single-frame** branches now call
   `omp_get_max_threads()` before `compute_ellipse()`. This does NOT help the
   case where the flag IS passed but with the value 1 (e.g. the Tcl side's own
   `resolve_job_count` resolving to 1, from a stale `parallel_jobs` setting or
   a cgroup-limited core count) — that is an explicit request for 1 thread and
   this binary honours it, by design: overriding an explicit
   `--asym-threads 1` would also silently break the batch-worker contract the
   day someone passes it there too. That failure mode, if it is what a slow
   deployment is hitting, needs a Tcl-side fix (`resolve_job_count` /
   `parallel_jobs`), not a C-side one.
   `--batch-asym-ellipse` is untouched — its branch never calls
   the new auto-detect helper, so `asym_threads` structurally stays at its `1`
   initializer regardless of what any caller does or doesn't pass (its own
   worker-process pool already parallelises across frames; auto-detecting
   there would make every one of N worker *processes* also spawn nproc
   *threads*, oversubscribing by nproc×).

Neither change touches `pa_de`, the Nelder-Mead simplex arithmetic, or its
iteration order, so no result can move by even one ULP from reassociation.

### Verified byte-identical

Built a reference binary straight from the pre-change source
(`git show HEAD:native/sos_triangle_fast.c`) and diffed its output against
the modified binary: `--asym-ellipse` and `--asym-ellipse-geo` on all 50 frames
of `traj50_bare_inputs` (498-517 spheres, 27,314-27,322 atoms per frame - real
sphere/atom-count diversity, not one frame repeated), single-frame vs
`--batch-asym-ellipse` (a 50-frame joblist in one process, and a 1-job list),
`--asym-threads 1` vs the new auto-detected thread count, and `OMP_NUM_THREADS`
overriding the auto-detect — all byte-identical (0 mismatches / 50 frames).
`vmdhole/tests/test_accel_parity.sh`'s CLI dispatch-chain probes (which
guard the exact `--asym-*` argument-parsing
class of bug fixed at `sos_triangle_fast.c:3072`) still pass against a binary
built from this source.

### Measured speedup (same machine, 16 cores)

Single frame (`frame_00000`, 508 slices, `--asym-threads 1` unless noted;
`/usr/bin/time`, 3 reps each):

| configuration | wall time/frame | vs. reported symptom |
|---|---|---|
| pre-change, `--asym-threads 1` (matches the user's ~3.6 s report) | 3.83-4.01 s | baseline |
| this change, `--asym-threads 1` (hoist only, same thread count) | 3.25-3.26 s | 1.19x |
| this change, no `--asym-threads` flag (new auto-detect) | 0.34-0.35 s | **11.2x** |
| pre-change, `--asym-threads 16` (threading was already there IF wired) | 0.43-0.44 s | 8.8x |

All 50 frames of `traj50_bare_inputs`, `--batch-asym-ellipse` in one process
(so per-worker cost, not wall - this is the number that matters for a batch
worker pool):

| binary | 50 frames | s/frame |
|---|---|---|
| pre-change | 192.4 s | 3.848 |
| this change (hoist only; batch stays single-threaded per worker by design) | 152.9 s | 3.058 |

1.26x, matching the single-frame hoist-only number above within noise.

10,000-frame extrapolation, all from the measurements above (rows marked
[50-frame mean] use the more robust `traj50_bare_inputs` figure; rows marked
[frame_00000] use the single-sample auto-threaded number, since that path was
only measured on one frame - see "Verified byte-identical" above for where the
50-frame coverage applies):

| path | time | note |
|---|---|---|
| pre-change, serial (matches the reported "10 hours") [50-frame mean] | 10.7 h | 3.848 s/frame × 10,000 |
| this change, serial (hoist only, still no threading engaged) [50-frame mean] | 8.5 h | 3.058 s/frame × 10,000 |
| this change, single-frame auto-threaded (this fix) [frame_00000 only] | 58 min | 0.345 s/frame × 10,000 |
| pre-change, `--batch-asym-ellipse`, 16 workers [50-frame mean] | 40 min | 3.848/16 × 10,000 s |
| this change, `--batch-asym-ellipse`, 16 workers [50-frame mean] | 32 min | 3.058/16 × 10,000 s |

**Diagnostic note, outside this file's scope:** 3.6 s/frame matches serial
execution almost exactly, which is the strongest evidence for *why* the field
report is slow. There are two distinct ways `_asym_ellipse_c`
(`vmdhole.tcl:26159`) can end up serial, and this fix only covers one of them:

1. `sos_triangle_has_feature asymthreads` returns false (a stale capability
   cache, or a deployed binary that predates the flag) → `--asym-threads` is
   never passed at all → **this change's auto-detect fixes it** (58 min for
   10,000 frames instead of 10.7 h).
2. `resolve_job_count` resolves to 1 (a `parallel_jobs` setting stuck at `1`,
   or a genuinely single-core-visible environment) → `--asym-threads 1` IS
   passed, explicitly → **not fixed by this change** (this binary honours an
   explicit request for 1 thread, by design - see change 2 above); only the
   hoist applies (10.7 h → 8.5 h). This needs a Tcl-side fix if it is what a
   slow deployment is hitting.

A one-line check in a live plugin session tells which branch applies:
`puts [::VMDHole::sos_triangle_has_feature asymthreads]` and
`puts [::VMDHole::resolve_job_count]`.

Also note: `~/hole2/exe/sos_triangle` (the binary path the plugin's Settings
dialog points at by default) was deliberately **not** rebuilt or replaced by
this change — file ownership for this task was `native/` only. Deploying
this fix requires rebuilding from this source and pointing the plugin at the
result (`native/build.sh`), a step outside this file's scope.

## Fortran patches (`connolly_patches/`)

These modify HOLE2's Fortran sources — `hole` and `sph_process` — and are
applied to a real HOLE2 checkout by `apply_patches.py` (they are not
standalone-buildable). Each is a modified copy of the upstream file of the same
base name; upstream's own modification history is preserved in the file headers.

| Patch | Upstream file | Modification |
|---|---|---|
| `coarea_fast.f` | `coarea.f` | 2D spatial grid over the active circles, built once per call and queried for every square in both the coarse pass and every refinement cycle, replacing a linear scan over all circles. |
| `hcapen_fast.f` | `hcapen.f` | Skips the full scan when the geometry is comfortably inside the safe bound, where the result is guaranteed to equal the scanned answer. |
| `holcal_par.f` | `holcal.f` | Two-pass CONNOLLY driver: pass 1 grows the whole centreline with `CONCAL` suppressed to discover every plane, so the expensive per-plane `CONCAL` work in pass 2 can be parallelised with OpenMP. Per-plane scratch files carry the PID and plane index so concurrent `hole` processes cannot collide; they are opened `STATUS='REPLACE'` (see note below). |
| `holeen_par.f` | `holeen.f` | Replaces the blanket `SAVE` — which made all locals static and therefore raced when `HOLEEN` was called from multiple OpenMP threads — with per-thread (`threadprivate`) caches. |
| `sphqpu_par.f` | `sphqpu.f` | Parallel dot-culling pass for `sph_process`: the cull runs in parallel into an accumulator mask, then the surviving dots are written serially in the original order. Falls back to the unchanged serial loop if the parallel path cannot be set up. |
| `machine_dep.g77` | same | RNG call-counter, so the seeded random sequence is reproducible across the parallel paths. |
| `Makefile` | same | OpenMP build flags. |

The parallel patches are ordered so that output does not depend on thread count
or scheduling: work is distributed for the *computation*, but results are
emitted in the original serial order. This is what makes the accelerated
binaries produce byte-identical *results* — the `.sph` file and every number —
which `verify.sh` checks. It does not extend to the whole `.out` text; the
next section says exactly where the line falls.

**Scope of that guarantee — what is and is not preserved.** Measured on 1BL8
protein-only at 1, 4 and 8 threads, for `conn` alone and `conn` + `capsule`, at
`shorto` 0 and 1:

| | preserved? |
|---|---|
| `.sph` file | **byte-identical to stock** in all four configurations |
| every profile number (radius, Requiv, area, centre) | **byte-identical** |
| `.out` text: dependence on thread count | **none** — 1, 4 and 8 threads give identical text |
| `.out` text: verbose per-plane CONNOLLY log | **not preserved** — see below |

The last row is a real limitation, not a rounding of the claim. The pass-1
prepass calls `CONCAL` with `SHORTO` hard-coded to 2, which silences every
`IF (SHORTO.LT.1)` / `IF (SHORTO.LT.2)` progress message inside `CONCAL` and
`COAREA` — "Connolly routine for this plane", "Have stored total of *n*
points", "Area calculation have found *n* spheres", the per-cycle "Connolly
area calc cycle" lines. Pass 2 does not re-run `CONCAL` for those planes; it
replays their `.sph` records from the per-plane scratch files. So that output is
never produced anywhere. Measured on `conn`, `shorto 1`: **215 lines that stock
prints are absent**, and only 2 lines are new. No number differs — every line
that appears in both is identical, and the `.sph` matches byte for byte.

Suppressing them is what makes the prepass parallel: emitted from worker
threads their interleaving would depend on scheduling. Preserving them would
mean buffering each plane's text in memory and replaying it in plane order —
the same restructuring already suggested for the scratch-file plumbing.

**One message is replayed**, because it is a result rather than progress: the
`initial point probe radius ... less than probe radius ... So using HOLE point
for calcs` notice, which reports that a plane fell back to the HOLE point.
`concal_par.f` records it (`CN_REC`/`CN_FLAG`/`CN_VAL`, threadprivate) instead
of printing when called from the prepass, and `holcal_par.f` replays the
recorded notices serially in ascending plane order. On `conn` + `capsule`,
`shorto 1`: 13 notices from stock and 13 from the patched build at 1, 4 and 8
threads, with identical values.

The recorded value matters. The notice reports `-NEWENG` from a `HOLEEN` call
made *inside* `CONCAL`, which is not equal to `holcal`'s own `STRRAD` for the
same plane — they are evaluated at different centres. Recomputing the condition
outside `CONCAL` selects a different set of planes (measured: 6 instead of 13),
so the value is carried out of the parallel region rather than re-derived.

> **Correction.** An earlier revision of this note claimed the whole `.out`
> text was byte-identical to stock, and described the ungated notice as a
> *theoretical* gap that "was not observed" and was left unpatched by design.
> All of that was wrong. It rested on a comparison made with a build whose
> patched objects had silently failed to compile, so both sides were stock. The
> `.out` was never byte-identical in any mode tested; the thread-ordering
> hazard was real and is now fixed; and the verbose-log gap above went
> unrecorded entirely.

> **`verify.sh` footgun:** Part E compares `$HOLE_EXE/hole` (default
> `~/hole2/exe`) against `$ACCEL`. Run `verify.sh` **before** installing to the
> live exe directory, or `$HOLE_EXE` is the accelerated build you just installed
> and Part E silently self-compares.

### Scratch-file open mode (portability, not a bug fix)

The per-plane scratch files in `holcal_par.f` are opened `STATUS='REPLACE'`
rather than `STATUS='UNKNOWN'`. Names include the PID and the plane index, and
the files are deleted once read, but an aborted run leaves them behind, so a
later run with a recycled PID could reopen one. `UNKNOWN` does not truncate on
reopen, which in principle would leave a stale tail for the pass-2 reader (which
reads to EOF) to pick up as real records.

**Measured: this is not reachable with gfortran**, which truncates a sequential
file at the last written record on `CLOSE` (a 45-byte stale file became 14 bytes
after writing one short record). The standard does not guarantee that for
`UNKNOWN`, so `REPLACE` is used to state the intent explicitly. The change was
verified byte-identical (`.sph` and `.out`) against the previously deployed
build on a seeded CONNOLLY run.
