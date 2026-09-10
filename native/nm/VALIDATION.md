# Validation status of the ported tools

Everything here is measured on the author's own Nav trajectory
(`traj_per/`, 18,677 protein atoms) against the HOLE run already stored in
`traj_per/hole_output_step5_assembly/`, using that run's own cards:
cpoint `-2.4354 -0.9136 3.8850`, cvect `0.0166 0.0289 0.9994`,
sample 0.25, endrad 15.0, `simple.rad`.

## nm_search — as accurate as HOLE is with itself; the search kernel 5x faster

Scope of the 5x: the search step alone, on the author's Nav trajectory, which
is not a benchmark fixture. It is not a benchmark-suite number. In the
plugin's 15-worker pool on the 50-frame benchmark trajectory the whole
calc deliverable is 0.98x the accelerated HOLE search
(`paper/benchmarks/results/endtoend.csv`, row `vmdpathfinder_nm`): the pool is
I/O-bound there and the search step is not what limits it.

| | HOLE | nm_search |
|---|---:|---:|
| search, mean of 10 frames | 70 ms | **13.6 ms** |
| bottleneck radius, 9 of 10 frames | | within 0.06 Å of HOLE at seed 1 |
| bottleneck radius, frame 60 | 1.200 Å at seed 1 | 0.767 Å |

Frame 60 looked like a defect and was investigated to the bottom. It is not
one. The chain of evidence, in order:

1. **Identical atom model.** HOLE's sphere centre at t = 6.75 evaluates to a
   clearance of 1.4315 Å in the port's own atoms and radii; the port's
   centre, 3.3 Å away in the same plane, evaluates to 0.767 Å. Both programs
   agree about every atom; they chose different points.
2. **Two basins, no wall.** Mapping the clearance field in that plane shows
   the port in a basin of 0.70 Å and HOLE in one of 1.37 Å, 3.25 Å apart,
   joined by a saddle of 0.51 Å. Two parallel passages; each is a true local
   maximum.
3. **Every deterministic fix trades one frame for another.** A fine-pass
   probe wide enough to reach the second basin (±3 Å) lifts frame 60 to
   1.201 Å - and hops frame 90 from 0.234 Å into a contiguous side channel
   at 1.115 Å. Letting the repair sweep keep a slice's better clearance does
   the same. An in-plane connectivity gate does not separate the cases (the
   frame 90 channel is connected, and 93 slices long), so persistence along
   t cannot either. A probe-free sequential march reproduces HOLE more
   closely on frames 10, 40 and 80 but still gives 0.767 on frame 60.
4. **HOLE is bimodal on exactly these frames.** Same cards, eight seeds:

   | frame | seeds 1-8 | port |
   |---|---|---:|
   | 60 | 1.200 0.690 0.690 1.200 1.200 1.200 1.200 0.560 | 0.767 |
   | 70 | 0.670 0.670 0.670 1.000 1.000 0.670 0.670 0.670 | 0.608 |
   | 80 | 1.390 x8 | 1.360 |
   | 90 | 0.230 0.230 0.230 0.230 0.230 0.220 0.220 0.230 | 0.234 |

   And by Monte Carlo step count at seed 1, frame 60 gives 0.75, 0.76,
   0.56, 1.20, 0.69 for 100, 300, 500, 1000, 3000 steps. The stored result
   of 1.200 Å for that frame is one outcome of a coin the plugin's fixed
   `raseed 1` happens to land the same way every time.

So the port is within HOLE's own run-to-run distribution on all ten frames,
and deterministic. Where the structure has one passage the two agree to a
few hundredths of an Ångström; where it has two, HOLE picks one at random
and the port picks the narrower one every time.

### What was tried and did not survive (measured, do not re-run blind)
- Fine-pass probe at ±0.75 Å: frame 60 worse (0.767 → 0.561). Removed.
- Fine-pass probe at the cascade's ±1.5/±3 Å, with or without an in-plane
  connectivity gate: fixes 60, hijacks 70/80/90. Off by default
  (`NM_FINE_GATE`, `NM_FINE_SCALE`).
- Repair sweep keeping the better clearance at narrow slices: same trade.
  Reverted to unconditional replace.
- Coarse cascade scale/rounds sweep (9 combinations), repair tolerance sweep
  (5 values), coarse sampling at 4x and 2x (worse: 0.689), grid cell size
  3-20 Å: no useful effect. Knobs kept as `NM_PROBE_SCALE`, `NM_PROBE_ROUNDS`,
  `NM_JUMPTOL`, `NM_CSAMP_MULT`, `NM_CELL`.
- Sequential fine pass with the cascade at every slice: worse on 60/70/80/90
  and 5-7x slower. Without the cascade (`NM_SEQ=1 NM_SEQ_PROBE=0`): closer
  to HOLE on frames 10/40/80, same on 60/70, worse on 90 (0.408); kept as an
  option, not the default.

### What DID matter
**Thread count, by a factor of 17.** See `nm_threads.h`. OpenMP's default is
`nproc`, which on this 8-core machine is 16 logical CPUs, and that setting
runs the search at 368 ms against 21 ms on 8. This was the entire gap
between "the port is 5x slower than HOLE" and "faster than HOLE".

**Two changes to the clearance loop, which gprof puts at 100% of the search
time.** Both leave every output byte-identical on all ten frames:

| | mean search |
|---|---:|
| after the thread fix | 19.0 ms |
| + squared-distance reject before the sqrt | 15.7 ms |
| + atoms stored in cell order, x/y/z/r interleaved | **13.6 ms** |

### Consequence for the dropdown
Offer both. Label Monte Carlo as HOLE's own stochastic search, the one the
published numbers came from, and Nelder-Mead as deterministic and 5x
faster. On a structure with parallel passages they can differ by up to the
width of HOLE's own seed spread; that is a property of the pore, and the
honest thing is to surface it, which the Monte Carlo knobs (seed, steps)
already let a user do.

## Vestibule runaway (`NM_MAXSTEP`) — flat-then-cliff mouths, fixed 2026-09-05

Reported against the plugin's own Pore Profile plot: with the Nelder-Mead
search, the bulk region at one or both ends of the trajectory sometimes
plotted as a flat run followed by a vertical jump, instead of the smooth
funnel HOLE's Monte Carlo search draws there.

Root cause, found by dumping the raw slice centres (not just the radius):
`nm2d_once`'s per-slice optimizer is a plain, unconstrained Nelder-Mead
simplex. `nm_cap` (endrad + margin) stops it once its VALUE is clearly
outside the pore, but a vestibule opening into bulk solvent is a nearly-flat
landscape over a wide area - every reflect/expand keeps finding a genuine,
if marginal, improvement in radius, so the value check never fires while the
position drifts arbitrarily far. This is the standard Nelder-Mead pathology
on a flat objective, not a wall crossed (`nm_connected()` would still pass -
the whole path is typically open). It is not new: it reproduces on the
already-validated 10-frame corpus (`fin_*.txt`), just not on the metric that
corpus checked (bottleneck radius, unaffected by a mouth-only artifact).

| frame | pre-fix max single-slice jump | new max (capped) |
|---:|---:|---:|
| 0  | 19.5 Å  (r 11.34 -> 18.52) | 6.0 Å |
| 10 | 24.0 Å  (r 10.78 -> 20.20) | 6.0 Å |
| 20 | 39.0 Å  (r  7.78 -> 22.40) | 6.0 Å |
| 80 | 38.0 Å  (r  9.41 -> 23.59) | 6.0 Å |
| 30, 50, 90 | 2.3-5.4 Å (already mild) | unchanged in shape |

Fix: `NM_MAXSTEP` (default 6 Å, matching the offset probe cascade's own
2-round combined reach) rejects any reflect/expand candidate farther than
that from where the call started, forcing the simplex to contract and
converge locally instead of running away. `nm2d_once` is called from nine
sites (coarse march, its probe cascade, the fine pass, its gated probe, the
repair sweep, capsule march); the cap applies uniformly since all nine are
documented as local refinement from an already-good seed, never the
mechanism meant to reach a distant basin (that is the explicit offset
probes, whose own max reach is comfortably under the cap).

Checked on the full 10-frame corpus after the fix:
- Bottleneck (min) radius: identical to 1e-5 on 9 frames, 0.00001 Å on the
  10th (frame 90) - the accuracy metric the corpus was built to check is
  untouched.
- Deep-pore rows (away from a correction): 63-97% byte-identical per frame;
  the rest differ by under 0.001 Å, propagated forward from the corrected
  mouth through the coarse march's carried-forward seed, not new noise (a
  cap set to a no-op value reproduces the pre-fix output exactly, confirming
  the code path is otherwise deterministic).
- `--conn`: unaffected by the fix directly (`concal`/`coarea` call
  `clearance()`, not `nm2d_once`), but its point cloud follows the corrected
  centreline, so a Connolly run on a fixed mouth gets ~50 more, smoothly
  spaced dots there instead of a few spanning a 20-40 Å teleport; Requiv on
  every shared centre is unchanged (< 0.01 Å).
- Visual: the plugin's own Pore Profile plot, same frame and CPOINT/CVECT
  that showed the bug, goes from a smooth rise into a flat plateau (r~7-8
  over ~7 Å) then a vertical jump to 18.5, to a smooth continuous rise with
  no plateau or cliff.
- `test_nm_engine` (37 checks), the 21-group main suite and the 14-test unit
  suite all still pass.

Smaller caps (2, 3, 4 Å) were tried and rejected: they still remove the
runaway, but perturb rows throughout the WHOLE profile, not just after a
correction (down to ~50-60% exact-match), because a cap that tight
occasionally clips legitimate multi-Å local convergence even deep in the
pore - the same convergence-noise sensitivity `nm2d_once`'s own header
comment already documents (a sub-0.01 Å seed difference once flipped an
entire track, frame 0 rms 0.027 -> 1.34 Å). 6 Å is the smallest value tested
that stays out of that regime.

## Connolly pass (`--conn`) — byte-identical to HOLE on the same centres

The plugin's Connolly surface is HOLE's `conn` card: per centreline sphere,
`concal.f` grows probe positions in the slice plane, `coarea.f` integrates the
union of their in-plane circles for Requiv, `addend.f` adds escape spheres at
both ends. The plugin already carries a bit-exact Tcl port of all three; the C
in `nm_holeout.h` is that port line for line.

Checks on frame 0 (18,677 atoms, `conn 1.15 0.5`, sample 0.25, endrad 15):

| comparison | result |
|---|---|
| C vs the Tcl port, six slices (-175, -139, -90, -89, 0, 20) on identical rounded centres | byte-identical `.sph` lines, all six |
| C on HOLE's own centres (read back from `hconn.sph`, 3 decimals) vs HOLE's dots | 108,303 vs 108,334 dots, 229 differ (0.2%); 273/315 slices identical in dot count, 281/315 identical Requiv; the rest differ by one probe at the 1.15 A boundary or 0.01-0.09 A Requiv - input rounding, since the Tcl port on the same rounded centres gives the C output exactly |
| plain loops vs accelerated (hash + block fill), HOLE centres and NM centres | byte-identical `.sph` and `.tsv` |
| plugin run on 1GRM, Monte Carlo centres + fast Connolly vs HOLE `conn` | 845 dots both, min Requiv 0.409 both (`test_nm_engine`) |

Time, same frame: HOLE spherical+conn 480 ms; port 336 ms single-thread and
55 ms on the physical cores (search 14 ms on top). The plain port was
1,127 ms single-thread - 53% in the area raster, 36% in the O(n^2)
duplicate test (298 M calls) - before the two exact accelerations.

One guard the Tcl port does not need: a slice whose centre clearance is already
past endrad+3 (NM's terminal slices can be; HOLE's search stops just past
endrad) would loop forever in `concal`, so it is treated as escaped.

## Thread policy under load — passive waiting, set by a one-time re-exec

The search runs one short OpenMP region per slice. With GOMP's default
spin-wait and 16 busy loops on the 16 logical cores, the 14 ms search took
474 ms at 8 threads (HOLE under the same load: 1.26 s). `OMP_WAIT_POLICY=
PASSIVE` gives 63 ms at 8 threads, 75-87 ms at 1-4. The runtime reads the
policy before `main()`, so `nm_set_wait_policy()` sets it and re-executes the
binary once (`/proc/self/exe`, then `argv[0]`); `NM_NO_REEXEC` stops the loop
and a caller-set `OMP_WAIT_POLICY` is respected. Idle cost: none measurable
(13.6-13.8 ms). In the plugin's pools (`njobs > 1`) the engine also gets
`OMP_NUM_THREADS=1`, one frame per core.

Plugin wall time per frame on the same structure, idle, display off / on:
Monte Carlo 170 / 157 ms, Nelder-Mead 96 / 89 ms - the rest of the pipeline
(PDB write, parsing, meshing) is shared, so the search's 5x shows up as about
2x per frame.

## Where the search speed shows up in the plugin (Nav, 197k-atom system, selection protein, 18,677 atoms)

Per-frame wall time through `run_analysis`, this 16-thread box:

| run | Monte Carlo | Nelder-Mead |
|---|---:|---:|
| 100 frames, spherical, 15 jobs | 29-32 ms | 30-32 ms |
| 30 frames, spherical, 4 jobs | 35 ms | 36 ms |
| 30 frames, spherical, 1 job | 109 ms | 69 ms |
| 100 frames, Connolly, 15 jobs | 167 ms (HOLE conn) | 52 ms (48 ms after the packed-record reader) |
| 100 frames, Connolly, MC centres + ported conn | 74 ms | - |

With several jobs the engines are hidden behind the pool and the wall time is
the plugin's serial work per frame, which a builtin-level profile puts at
~20 ms in `open "|sh run.sh"` alone (808 opens, 2.1 s of a 2.9 s run): forking
the VMD process, whose address space holds the whole trajectory. Everything
else per frame (coordinate record 1.6 ms, parse 1.7 ms, control file, axis)
is ~4 ms. So on a many-core machine a spherical trajectory run is fork-bound
whichever engine searches, and the search's 5x only appears with few jobs or
under Connolly. One persistent shell per pool slot instead of one fork per
frame (now in the run loop, `_pw_*` procs) halves it: 100 Nav frames at 15
jobs 30 -> 15 ms per frame for both engines, Nelder-Mead Connolly 48 -> 41 ms,
every frame's `.sph` and `.tsv` byte-identical to the per-frame launcher, and
killing all workers mid-run loses only the jobs in flight.

The engine reads the packed coordinate record (`input_frame.vhb`) directly,
rounding through `%.3f` so it matches the PDB path byte for byte
(`test_nm_engine` checks it); before that it cost an extra 13 ms `writepdb`
per frame.

## mesh_csg — wired as the spherical mesher

Vertex distance to the exact clipped sphere union, Nav frame 0 (318 pore
spheres, 185 clip spheres), every vertex:

| grid | triangles | mean | p95 | area vs 0.5 Å | build |
|---|---:|---:|---:|---:|---:|
| 1.0 Å uniform | 8,357 | 0.0053 Å | 0.035 Å | -2.0% | 16 ms |
| 0.75 Å uniform | 15,087 | 0.0032 Å | 0.018 Å | -0.9% | 31 ms |
| 1.0 Å + 0.5 Å neck | 20,845 | 0.0026 Å | 0.013 Å | -1.7% | 30 ms |
| **1.4 Å + 0.7 Å neck (shipped)** | **12,238** | **0.0044 Å** | **0.023 Å** | **-2.9%** | **17 ms** |
| 0.9 Å + 0.45 Å neck | 23,281 | 0.0029 Å | 0.015 Å | -1.4% | 34 ms |
| 0.5 Å uniform | 34,199 | 0.0015 Å | 0.007 Å | 0 | 86 ms |

Solving the edge crossing against the governing sphere instead of
interpolating the sampled field halves the mean error at every spacing
(1.0 Å: 0.0100 -> 0.0053; 1.0/0.5: 0.0054 -> 0.0026) for no extra time.

The maxima sit where a clip sphere cuts a pore sphere (the field is a
CSG max, not a true distance there). Playback on the 10-frame Nav
trajectory, plugin cost per newly visited frame, Xvfb:

| | sos_triangle (stride-4 draft) | mesh_csg (1.0/0.5 Å refined, full mesh) |
|---|---:|---:|
| first visit | 185 ms (86 sph_process + 78 sos_triangle) | 130 ms (16 mesh + 110 parse/draw of 21k triangles) |
| revisit (cached file) | 7 ms (1.7k triangles) | 77 ms (21k triangles) |
| settle | 25 ms (same 7k-triangle mesh) | 211 ms (34k triangles at 0.5 Å) |

One mesh now serves both playback and the settled view, so the surface
cannot change shape as it settles. Two coarse-then-fine schemes were tried
first and both were rejected on that alone: a uniform 1.0 Å playback mesh
(the neck visibly changed shape), then a 1.0/0.5 playback mesh against a
0.5 Å settle (the wide regions changed). Sizing the single mesh is then a
budget question - drawing costs about 2.9 us per triangle - and 1.0/0.5 is
the most accurate mesh that fits: 84 ms per newly visited frame on a
200k-atom, 10-frame trajectory (20 ms mesh, 60 ms draw), against 118 ms for
0.8/0.4 and 98 ms for 0.9/0.45. A smaller refinement radius is a false
economy: R=2 Å drops 1.0/0.5 to 8.4k triangles but takes the neck error
from 0.0019 to 0.0069 Å, because the neck's own spheres are 2-3 Å.
Two write-side fixes came out of the same work: the plot writer's printf
was half the mesh time and is now a fixed-point formatter (io 38 -> 6 ms),
and the records address a molecule directly so the plugin can source them.

The shipped default was then set from the granularity `sos_triangle`
actually produces rather than from the error columns. At dot density 15 its
median triangle edge is 0.86 Å; 1.4/0.7 gives 0.74 Å, so the two surfaces
look alike, where the earlier 1.0/0.5 default was 0.54 Å - three and a half
times more triangles than the pipeline it replaced. Over that whole range
the mean vertex error only moves from 0.0026 to 0.0044 Å, far below
anything the geometry means, so triangle count, not error, is the honest
knob. Nav playback at the shipped default: 51 ms per newly visited frame
(12 ms mesh, 35 ms draw, 11.4k triangles).

## Tunnels

A tunnel .sph is a plain sphere union with no clip records, so the same
mesher serves it. Meshing goes through the session's `--serve` child rather
than the job pool: process creation from a loaded VMD costs ~16 ms, more
than meshing one tunnel takes. KcsA, four routes, mesh plus draw:

| | sos_triangle | mesh_csg |
|---|---:|---:|
| mesh + draw | 217 ms | 59 ms |
| mesh alone | 189 ms | 30 ms |
| primitives drawn | 8,352 | 8,836 |

A property-coloured route still uses `sos_triangle`, for the reason below.
The two meshers write different filenames, so switching rebuilds rather
than serving whichever mesh is newer than the .sph.

## Connolly

HOLE's Connolly output is a union of spheres too - its `S-999` records carry
a radius (the 1.15 A probe up to the pore spheres), which is why sph_process
expands each one into a patch of dots. Sampled vertices of both meshes sit on
the same analytic surface: sos_triangle 0.0000 A from it, mesh_csg 0.0042 A.
So Connolly needs no new algorithm.

The grid is applied uniformly there: a molecular surface is at probe scale
everywhere, so refining by sphere radius would refine all of it. Nav, one
Connolly frame, plugin cost per newly visited frame:

| | sos_triangle | mesh_csg 1.4 A |
|---|---:|---:|
| plain colouring, first visit | 79 ms | 63 ms |
| lobe colouring, first visit | 94 ms | 64 ms |
| lobe colouring, settle | 280 ms | 53 ms |

The settle is where it pays: sos_triangle draws a reduced mesh while playing
and the full one at rest, and the mesher draws one mesh for both.

Emission is per z-slab rather than per thread, so a given input meshes
byte-identically whatever the schedule and the thread count do. That is not
cosmetic: `conn_lobes split` memoises the first triangle to land in each
classification cell, so a mesh whose triangle ORDER moved between runs
coloured a few lobe triangles differently each time. It cost the
native-vs-Tcl byte-identity check in `test_conn_lobes_engine` before the
per-slab buffers went in.

HOLE marks a CUTTER by writing `LAST-REC-END` after the record, and
sph_process drops those from its surface passes (sos_triangle then culls the
triangles that touch them). mesh_csg reads that marker: in a spherical .sph
all 185 marked records also carry beta 0, so the older beta-only rule agreed
there, but a Connolly render file marks 1,874 records that keep beta 999.99.
Meshing those as surface is what put blobs on the ends of the lateral
openings; honouring the marker takes that frame from 42,852 to 14,038
triangles and matches the sos_triangle surface again.

Two caveats. The mesh is written with `--draw`, because `conn_lobes split`
finds triangles by the literal "draw trinorm" - the region split itself is
unaffected, since it assigns each triangle by its nearest classified dot, so
HOLE's own dots keep doing the classifying. And the marching-cubes surface
carries 4-5% less area than sos_triangle's, which does NOT converge away with
resolution (-4.7% at 0.5 A, -4.2% at 0.3 A with 138k triangles): the surface
has sharp cusp lines where probe spheres meet and marching cubes rounds them.
Side by side the two look equivalent, the marching-cubes one slightly smoother.

## The profile stops at ENDRAD, like HOLE's

The port used to keep the slice that crossed ENDRAD in the profile as well as
in the .sph. HOLE keeps it only in the .sph (it becomes a mouth clip sphere).
Three such rows, up to 17.18 A against ENDRAD 15, moved the r^2-weighted
Volume readout by 20% on a pore whose bottleneck agreed to 0.003%. Dropped
from the profile rows now; same frame, trapezoid of pi r^2 over the trace:
HOLE 9,543 A^3, port 9,508 A^3.

## Sourced meshes and the material

The addressed-form file began with `graphics $mol delete all`. Sourcing it
ran that line AFTER the plugin had set the material, and VMD's delete-all
discards the material entries, so every sourced surface drew in the default
material (the parsed path skips that line for exactly this reason). The
addressed form no longer carries it; the plugin clears the molecule itself.
Verified: the first primitives on a sourced surface are now `materials` and
`material`. The cache recipe token moved (`csg` -> `csg2`) so meshes written
before the change rebuild instead of being served.

## Property colouring: profiled, and the mesher is not the lever

Measured per newly visited frame, 9 Nav frames, property colouring on:

| | ms |
|---|---:|
| pore-lining sidecar (`write_hydro3d_atoms_sidecar`) | 33 |
| recolour (`sos_triangle --recolor`, one process per frame) | 36 |
| draw | 16 |
| **mesh build** | **absent - the base mesh is cached** |

Meshing is not in the picture, so a faster mesher cannot help, and a denser
one hurts: the recolour and the draw both scale with triangle count. That is
the whole reason the wiring below was reverted.

The two real levers are the 33 ms sidecar, which is a Tcl loop over the .sph
plus atomselect work, and the ~16 ms of the recolour that is process creation
from a loaded VMD (sos_triangle already has `--batch-hydro3d-recolor` for the
multi-frame case; the interactive path has no equivalent).

## The Connolly property failure was the prebuild, not a polygon limit

"sos_triangle produced no surface (likely its polygon limit...)" on a
property-coloured Connolly run of a 200k-atom system was not sos_triangle
at all. `prebuild_surfaces_parallel` skips a large Connolly cloud on purpose
(the lazy path trims and reduces it first), but its registration loop then
marked every frame with no mesh as FAILED with that message, and
`load_surface_for_frame` refuses a failed asset without retrying. run_analysis
calls the prebuild whenever property colouring is on, so the frames were dead
before they were ever drawn - with either mesher. Skipped frames are now left
for the lazy path. Nav, Connolly, property colouring, per newly visited
frame with the mesher: 73 ms first visit, 20 ms revisit, 14,038 coloured
triangles.

## Property colouring: adopted after all, for a different reason

The first attempt fed the mesher's mesh to sos_triangle's recolour and lost
on time (149 ms against 93). It came back for two reasons that outweigh
that: the recolour has no polygon limit on the mesher's mesh, which is what
let a dense Connolly run colour at all, and the default grid is now coarser
(1.4/0.7), so the recolour has fewer triangles to colour. Measured, Nav,
spherical, property colouring, per newly visited frame: sos_triangle base
89 ms, mesher base **22 ms**; settle 25 -> 10 ms.

## Property colouring: possible, measured, not adopted (first attempt, kept for the record)

`--draw` writes HOLE's own records, which is what sos_triangle's
per-triangle recolour modes read, so a marching-cubes mesh can be coloured
by a property. It was wired up and verified correct on the Nav trajectory:
7 distinct colours over 2,172 colour runs on the mesher's mesh against 621
on sos_triangle's, same scheme, same range.

It was then removed, because it is slower:

| property colouring, per new frame | sos_triangle | mesh_csg |
|---|---:|---:|
| total | 93 ms | 149 ms |
| triangles coloured and drawn | 6,018 | 12,453 |

The recolour costs about 9 us per triangle, several times what drawing
costs, and marching cubes emits roughly twice the triangles sos_triangle
does at the same median edge length (it produces a long tail of slivers).
Meshing was never the expensive half of a property frame, so making it
faster cannot pay for a denser mesh. Colouring inside the mesher, which
would skip the recolour pass entirely, is the only version of this that
could win.

Not done at that point: mesh_csg had no capsule or Connolly path (their
`.sph` files are not sphere unions).

## One shipped binary

`sos_triangle_fast` now carries nm_search, mesh_csg and conn_lobes
(`-DVMDPATHFINDER_MULTICALL`; each tool's `main` becomes `<tool>_main` and the
dispatcher routes `--nm-search`, `--mesh` and `--conn-lobes`). One trap:
nm_search re-execs itself once to pin `OMP_WAIT_POLICY=PASSIVE`, with the
argv it was handed - the shifted one, minus the subcommand - so the re-exec
ran the triangulator. The dispatcher does that re-exec itself with the full
argv. With the three standalone files moved out of the exe directory, the
mesher, search and classifier test groups all pass through the one file.

## Property colouring in-process, and what a property frame costs now

| per newly visited frame, Nav, 9 frames | sos recolour | mesher recolour |
|---|---:|---:|
| spherical, property | 149 ms | 103 ms |
| Connolly, property | 468 ms | 325 ms |
| recolour alone | 44 / 141 ms | 8 / 11 ms |
| lining sidecar, Connolly | 127 ms | 33 ms |

The recolour is now a persistent-process call and the sidecar no longer
parses the .sph in Tcl (the mesher reports the centreline extent). What
remains of a property frame is the atom-selection loop that writes the
sidecar (~30 ms) and the mesh build itself.

## Dots display

| per newly visited frame | sph_process + sos_triangle --points | mesher --dots |
|---|---:|---:|
| spherical | 19 ms | 3 ms |
| Connolly | 68 ms | 26 ms |

## One surface pipeline

Every mesh the plugin draws now goes through three procs in `vmdpathfinder.tcl`:
`surface_plot_name` (the one naming rule), `surface_mesh` (build one plot
from one .sph: the mesher when it can serve the run, the legacy
sph_process + sos_triangle pair otherwise and as the fallback, with the
sos dot-budget retry inside it) and `surface_mesh_cmd` (the same build as a
shell command for the job pools). The pore surface, every Connolly region,
every tunnel, the mean tubes, the dots display and the base a property
colouring starts from all call them; a smoke check now asserts that no
other proc runs sph_process or sos_triangle for a mesh. Before this the
same work was spread over thirty call sites, five mechanisms and six naming
rules. Tunnel meshes changed name in passing (`tunnel_NN.vmd_plot`, tagged
under the mesher); old `tunnel_NN.plot` files are simply rebuilt once.

## Capsule

Upstream HOLE writes a capsule run's slices as QC1/QC2 record pairs (the two
cap centres, radius in occupancy, beta 0.00) and `sph_process` has its own
capsule pass (`sphqpc.f`). That pass also meshes the escaped search records
(`-888`, radii to 26 A), so on the frozen 1GRM capsule fixture the legacy
surface is 7758 triangles spanning z -54..49 A for a pore whose kept cap
centres lie in z 2.8..5.1 A. The plugin used to draw a Tcl stadium tube
instead. Now the mesher reads a QC pair as one capsule primitive (distance to
the segment minus R; a sphere is the degenerate case) and takes `--axis
CPOINT CVECT ENDRAD`: a slice is kept only when both cap centres lie within
ENDRAD of the axis and its equal-area radius does not exceed ENDRAD, the rule
HOLE's own profile applies and the rule the Tcl tube used. The run's
escaped records are dropped, not used as cutters: as cutters they removed
60% of the tube.

| fixture, 1.0/0.5 A | slices kept | triangles | z extent | time |
|---|---|---|---|---|
| mesher, no axis | 67 of 67 | 16008 | -4.4..11.0 A (escaped slices 25 A long) | 14 ms |
| mesher, `--axis 0 0 0 0 0 1 8` | 32 of 67 (= the Tcl rule) | 1940 | 1.8..6.1 A | 9 ms |

Under the sos mesher the same kept slices go to `sph_process`, whose capsule
pass (`sphqpc.f`) then draws the tube instead of the blobs: 406 triangles,
z 1.9..6.0 A at dot density 4. The pure-Tcl fallback carries a port of that
pass (sphqpc.f, hocapr.f, hocapd.f, and SPHCHC's pass headers, including
HOLE's quirk that under `-sos -colour` the high-radius band goes out under
the colour -1 end-cap header); its `.sos` is byte-identical to the binary's
on the kept fixture, plain and coloured, and the spherical output is
unchanged.

Live capsule run on 1GRM (virtual display): 5172 triangles, wireframe 7753
lines, dots 2586 points, centreline 106 segments (two tracks), property
colouring through the same recolour path. The spherical, Connolly, lobe
and search outputs are byte-identical before and after (13 comparisons).

## Surface smoothing

Two local averages of the surface, one per mesher, behind one plugin entry
(`surface_mesh ... with`). Marching mesher, `--with`: the window frames' fields
on one grid, averaged, marched with linear edge crossings (the analytic
sphere crossing has no meaning on a mean field); caps and colour bands from
the centre frame. sos_triangle `--sos-smooth`: each centre dot moves to the
mean of itself and its nearest same-facing dot (normals agreeing, within
rho = 2 A) in every window cloud, the normal is the renormalised mean, header
and centreline records are copied; the pure-Tcl port is byte-identical.

Nav frames 0-2 (spherical, MC), centre = frame 1:

| path | single frame | window +-1 |
|---|---|---|
| mesher 1.4/0.7 A | 13146 tris, 15 ms | 12149 tris, 30 ms; mouths pulled in (z 47.6 -> 42.3 A, y 26.1 -> 19.8 A), body unchanged |
| mesher, the same frame as its own window x3 | 13146 tris | 13158 tris, 83% of vertices coincide within 0.1 A |
| sos pair, dot density 8 | 1742 dots, 2241 tris | 1742 dots moved 0.35 A mean / 1.22 A max (597 unmoved), 2473 tris |

The no-window mesher path is byte-identical to the previous binary.

## The pore-lining sidecar

With the mesher colouring, the sidecar is one C-speed `writepdb` of the
boxed selection with each residue's value in the B-factor column and a 1.0
in occupancy on the CA atoms (the molecule's own beta/occupancy are put
back afterwards); the mesher's `--atoms` loader takes either layout. Same
frame, same box: 76 ms for the Tcl per-atom loop, 17 ms for the PDB write;
13,146 of 13,146 triangles coloured identically from either. The sos batch
recolour prebuilds are skipped under the mesher, since the per-frame
recolour they were optimising is now 8-11 ms in-process.
