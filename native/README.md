# `native/`: the compiled engines and accelerators

Everything in here is native (C / patched Fortran) code the plugin can use:
the fast `sos_triangle` replacement, the MOLE 2 tunnel-search engine
(`mole_tunnel_engine`), the hydration accelerator (`hydro_project`), and the
optional HOLE-side Fortran patches.

`sos_triangle_fast.c` is a **modified copy** of HOLE 2.x's `sos_triangle.c`
(from <https://github.com/osmart/hole2>). It is a **drop-in replacement** —
same CLI, same input and output formats, **byte-for-byte identical** surfaces,
just much faster (`./verify.sh` reproduces that claim; measurements are in the
paper). It also adds an opt-in hydrophobicity-colouring capability the plugin
uses automatically when present; see [`CHANGES.md`](CHANGES.md) for the full
statement of modifications.

## Quick start

**Install everything, in one step:**

```sh
./build-vmdhole-optimized.sh                                       # build into build/
HOLE_INSTALL_DIR="$HOME/hole2/exe" ./build-vmdhole-optimized.sh    # build and install
```

It builds every binary the plugin uses: HOLE, `sph_process`, the fast
`sos_triangle`, and `mole_tunnel_engine` (cloning HOLE2 if you don't pass a
source path). It checks the toolchain first and backs up whatever it replaces.
Without `HOLE_INSTALL_DIR`, point VMDHole at `build/` via
**Settings… → Engines**.

## What's here

| path | purpose |
|---|---|
| `sos_triangle_fast.c` | the optimised, drop-in `sos_triangle` source |
| `voronoi/` | exact-arithmetic Delaunay/Voronoi predicates shared by `sos_triangle_fast` and the MOLE engine, with their unit tests |
| `xalloc.h` | checked-allocation helper shared by the C above |
| `mole/` | the MOLE 2 tunnel-search port and its per-unit test tools |
| `hydration/` | the hydration accelerator (`hydro_project.c`) with its reference implementation and tests |
| `connolly_patches/` | the optional HOLE-side Fortran patches (`apply_patches.py` + the patched units) |
| `build.sh` | builds the three self-contained C tools — what CI runs on a fresh checkout, no HOLE tree needed |
| `build-vmdhole-optimized.sh` | the full install script: patches and rebuilds HOLE + `sph_process` from a HOLE source tree and builds the C tools |
| `verify.sh` | rebuild upstream + fast `sos_triangle` and confirm byte-identical output |
| `CHANGES.md` | the statement of modifications relative to upstream |
| `LICENSE`, `NOTICE` | Apache-2.0 license and attribution for this derivative work |

Compiled binaries land at this directory's top level and are never committed.

## License

This directory is a derivative of Apache-2.0 software and is distributed under
the **Apache License, Version 2.0** (`LICENSE`), with attribution in `NOTICE`
and a modification notice at the top of `sos_triangle_fast.c`. This is
separate from the VMDHole plugin's own license (see the repository root).
HOLE 2.x © G. M. P. Coates, O. S. Smart, SmartSci Limited.
