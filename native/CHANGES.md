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
| `hcapen_fast.f` | `hcapen.f` | Skips the full scan when the geometry is inside the safe bound, where the result is guaranteed equal. |
| `holcal_par.f` | `holcal.f` | Two-pass CONNOLLY driver so per-plane work can run under OpenMP; per-plane scratch files carry PID + plane index, opened `STATUS='REPLACE'`. |
| `holeen_par.f` | `holeen.f` | Blanket `SAVE` (which raced under OpenMP) replaced with per-thread caches. |
| `sphqpu_par.f` | `sphqpu.f` | Parallel dot-culling pass; results written serially in the original order, with a serial fallback. |
| `machine_dep.g77` | same | RNG call-counter so the seeded sequence is reproducible across the parallel paths. |
| `Makefile` | same | OpenMP build flags. |

The parallel patches distribute the computation but emit results in the
original serial order, so output does not depend on thread count or
scheduling.

> Note: `verify.sh` Part E compares `$HOLE_EXE/hole` (default `~/hole2/exe`)
> against `$ACCEL`. Run it **before** installing to the live exe directory,
> or Part E silently compares the accelerated build against itself.
