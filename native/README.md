# `native/`: the compiled engines and accelerators

Everything in here is native (C / patched Fortran) code the plugin can use:
the fast `sos_triangle` replacement, the MOLE 2 tunnel-search engine
(`mole_tunnel_engine`), the hydration accelerator (`hydro_project`), and the
optional HOLE-side Fortran patches. `sh build.sh` builds the self-contained C
tools; `build-vmdhole-optimized.sh` additionally patches a HOLE source tree.

`sos_triangle_fast.c` is a **modified copy** of HOLE 2.x's `sos_triangle.c`
(from <https://github.com/osmart/hole2>). The upstream file is **not** modified.
The command-line interface, the stdin `.sos` input, and the stdout `.vmd_plot`
output are unchanged, so it is a **drop-in replacement** that produces
**byte-for-byte identical** surfaces, just much faster.

> `sos_triangle` is C, not Fortran (only `hole`/`sph_process` are Fortran).

It also adds one **opt-in** capability the upstream lacks: colouring the pore
surface by residue **hydrophobicity** in the same compiled pass (used
automatically by VMDHole when present; see [`CHANGES.md`](CHANGES.md)). This is
inert unless the `--hydro-atoms` flag is given, so the byte-identical drop-in
guarantee for normal use is unchanged.

## Quick start

**Install everything, in one step:**

```sh
./build-vmdhole-optimized.sh                                       # build into build/
HOLE_INSTALL_DIR="$HOME/hole2/exe" ./build-vmdhole-optimized.sh    # build and install
```

It builds every binary the plugin uses: HOLE, `sph_process`, the fast
`sos_triangle`, and `mole_tunnel_engine`. It checks the toolchain first and backs
up whatever it replaces. Without `HOLE_INSTALL_DIR`, point VMDHole at `build/`
via **Settings… → Engines**.

**Maximum speedup** (all HOLE binaries `-O2` *and* the fast `sos_triangle`,
collected into `build/`):

```sh
./build-vmdhole-optimized.sh [/path/to/hole2/src]   # clones HOLE2 if no path given
# then point VMDHole's hole / sph_process / sos_triangle paths at native/build/
```

## What's here

| path | purpose |
|---|---|
| `sos_triangle_fast.c` | the optimised, drop-in `sos_triangle` source (top level: it is the primary deliverable) |
| `voronoi/` | exact-arithmetic Delaunay/Voronoi predicates shared by `sos_triangle_fast` and the MOLE engine, with their unit tests |
| `xalloc.h` | checked-allocation helper shared by the C above |
| `mole/` | the MOLE 2 tunnel-search port: engine sources (`mole_main.c`, `mole_tunnel.c`, `mole_lining.c`, `mole_complex.c`, `mole_dh.c`, `mole_rng.c` + headers) and every per-unit test/dump tool the test suite compiles from source |
| `hydration/` | the hydration accelerator: `hydro_project.c` plus its reference implementation, generator and tests |
| `connolly_patches/` | the optional HOLE-side Fortran patches (`apply_patches.py` + the patched units) and their cache test |
| `build.sh` | builds the three self-contained C tools (`sos_triangle_fast`, `mole_tunnel_engine`, `hydro_project`) into this directory - what CI runs on a fresh checkout, no HOLE tree needed |
| `build-vmdhole-optimized.sh` | the full install script: patches and rebuilds HOLE + `sph_process` from a HOLE source tree AND builds the C tools; with `HOLE_INSTALL_DIR=...` it copies everything into your exe dir. (Rebuilding HOLE at `-O2` alone measures **2.04× (circular)** / **2.91× (CONNOLLY)** on 1BL8 with byte-identical output - HOLE's own Makefile sets no optimisation flag.) |
| `verify.sh` | rebuild upstream + fast `sos_triangle` and confirm byte-identical output (fetches the pinned upstream source if no local copy exists) |
| `CHANGES.md` | the six `sos_triangle` changes, complexity table, correctness argument |
| `LICENSE`, `NOTICE` | Apache-2.0 license and attribution for this derivative work |

Compiled binaries land at this directory's top level and are never committed.

## Results in one line

Identical output to upstream, much faster: **`sos_triangle` itself is ~9.4×
faster at the plugin's default dot density (15)**, growing to **~75×** on large
surfaces (the change is algorithmic — two O(N²) bottlenecks become
near-linear), and it raises the polygon cap so large pores that previously
aborted now complete. `hole` and `sph_process` are **also** accelerated
(Fortran patches, see [`CHANGES.md`](CHANGES.md)); the measured
**whole-pipeline** speedup against stock HOLE 2.2 at the same `-O2` is
**11.8× (circular)** / **6.9× (CONNOLLY)**. Full measurements, caveats, and
the data of record: [`CHANGES.md`](CHANGES.md) and
`../paper/benchmarks/results/`.

```sh
./verify.sh   # reproduce the "byte-for-byte identical" claim
```

## License

This directory is a derivative of Apache-2.0 software and is distributed under
the **Apache License, Version 2.0** (`LICENSE`), with attribution in `NOTICE`
and a per-Section-4(b) modification notice at the top of `sos_triangle_fast.c`.
This is separate from the VMDHole plugin's own license (see the repository
root). HOLE 2.x © G. M. P. Coates, O. S. Smart, SmartSci Limited.
