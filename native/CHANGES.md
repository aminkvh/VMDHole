# Changes relative to upstream HOLE2

The Apache-2.0 statement of modifications for `native/`: the C program
`sos_triangle_fast.c`, and the Fortran patches in `connolly_patches/` (applied
to a HOLE2 checkout by `apply_patches.py`). All modifications preserve output:
the accelerated binaries are byte-for-byte identical to upstream on the same
input — `./verify.sh` reproduces this claim. Rationale and detail are in the
source comments; measurements are in the paper (`paper/benchmarks/`).

## `sos_triangle_fast.c`

A modified copy of `sos_triangle.c` from HOLE2
(https://github.com/osmart/hole2, Apache-2.0). Drop-in replacement — same
CLI, same `.sos` input, same `.vmd_plot` output.

Performance/capacity changes (algorithm and output unchanged):

- `cull_coords`: spatial hash grid replaces the exhaustive pairwise duplicate
  test.
- `check_point`: O(1) used-vertex flag replaces a rescan of every triangle.
- `neighbour`: edge-invariant terms hoisted out of the candidate loop, and a
  uniform spatial grid restricts the search to the only dots that can pass
  the existing 3·|base| cutoff (ties still break by smallest dot index; a
  mis-estimated cell size falls back to the original full scan).
- `destroy`: hash on the endpoint pair replaces a linear scan of the edge
  list, with bucket chains in creation order so the first match is unchanged.
- `MAX_COORD` raised 30000 → 200000, so large/high-density surfaces that
  previously aborted now complete.

Opt-in additive features (inert without their flags; a normal `-s`/`-v`
invocation is still byte-for-byte upstream):

- `--hydro-atoms` / `--hydro-sph` / `--hydro-scheme kd|ww` /
  `--hydro-signed`: colour the surface by nearest-lining-residue
  hydrophobicity in the same compiled pass, bit-identical to the plugin's Tcl
  path. `--hole-features` is the capability probe.
- `--points`: emit unique vertices instead of triangles (the plugin's dots
  display).
- `--recolor` + `--hydro-values`, and `--batch-recolor FILE` /
  `--batch FILE` / `--batch-asym-ellipse FILE`: recolour or process many
  frames in one process.
- `--asym-ellipse` / `--asym-ellipse-geo` / `--asym-threads N`: the
  PoreAnalyser-style asymmetry fit; single-frame runs auto-detect the thread
  count when `--asym-threads` is not given.

Robustness fixes (no effect on valid input): the `calc_tri` recursion runs on
a dedicated large-stack thread with a depth ceiling instead of crashing on
very large meshes; a bare `-X` with no operand is reported instead of reading
past `argv`; `read_cord` checks its `sscanf` field count instead of silently
inheriting stale values from malformed records.

## Fortran patches (`connolly_patches/`)

Each is a modified copy of the upstream file of the same base name; upstream's
own modification history is preserved in the file headers.

| Patch | Upstream file | Modification |
|---|---|---|
| `coarea_fast.f` | `coarea.f` | 2D spatial grid over the active circles replaces a linear scan. |
| `hcapen_fast.f` | `hcapen.f` | The CAPSULE energy routine. A uniform grid over the atoms (built once per dataset) replaces the per-call scan of every atom: the cells around the query segment are gathered, a bound proves nothing outside them can rank in the two best, and the original loop is then replayed in original atom order over just those atoms - same expressions, same floating-point path, same result. Falls back to the full scan whenever the bound cannot be established. |
| `holcal_par.f` | `holcal.f` | Two-pass CONNOLLY driver so per-plane work can run under OpenMP; per-plane scratch files carry PID + plane index, opened `STATUS='REPLACE'`. |
| `holeen_par.f` | `holeen.f` | Blanket `SAVE` (which raced under OpenMP) replaced with per-thread caches. |
| `sphqpu_par.f` | `sphqpu.f` | Parallel dot-culling pass; results written serially in the original order, with a serial fallback. |
| `machine_dep.g77` | same | RNG call-counter so the seeded sequence is reproducible across the parallel paths. |
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

> The full `.out` text is not byte-identical to stock in any mode tested — the
> verbose-log gap above is real, not theoretical. A prior comparison that
> reported full text identity ran against a build whose patched objects had
> silently failed to compile, so both sides were stock.

> Note: `verify.sh` Part E compares `$HOLE_EXE/hole` (default `~/hole2/exe`)
> against `$ACCEL`. Run it **before** installing to the live exe directory,
> or Part E silently compares the accelerated build against itself.
