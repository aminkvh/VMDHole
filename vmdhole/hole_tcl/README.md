# hole.tcl — HOLE 2 in pure Tcl (standalone)

A from-source replication of HOLE 2 with no compiled dependency. Standalone so
it can be validated on its own; the same engine is inlined into `vmdhole.tcl`
as the plugin's pure-Tcl fallback.

Run it:

```sh
tclsh hole.tcl -pdb 1GRM.pdb -rad simple.rad -cvect "0 0 1" \
               -sample 0.5 -endrad 8.0 -csv profile.csv -sph out.sph
```

Defaults match the real binary's own banner: `MCSTEP 1000`, `MCLEN 0.1`,
`kT 0.1`, `sample 0.25`, `endrad 22.0`.

It covers the spherical, CAPSULE and CONNOLLY calculations and the 3D surface
pipeline, and refuses anything it cannot reproduce faithfully rather than
approximating it. It is much slower than the real binaries; the plugin only
uses it when they are unavailable.

## Validation (the HOLCAL trace)

Output is validated against seed-pinned reference builds of the real
binaries: `.sph` files byte-identical, profiles bit-exact, established by
step-by-step full-precision tracing of the annealing loop (the HOLCAL trace).
Per-file derivation and fidelity notes are in each source file's own header,
and `tests/` holds the regression guards.

## Licensing

HOLE 2 is Apache License 2.0 (`hole2/LICENSE`), so both this derivation and
redistribution of the original binaries are permitted.

## Reference binaries

The fidelity tests compare against `hole`, `sph_process`, `sos_triangle`
built from `native/stock_build/hole2/src` (upstream `osmart/hole2` commit
`a8eaf6121ba66625446933f4acd7d6aa336dbb47`), pristine — no local patches:

```sh
cd native/stock_build/hole2/src
source ../source.apache      # sets HoleVersion=2.3.1 etc.
make ../exe/hole ../exe/sph_process ../exe/sos_triangle
```

then copy the three binaries into `vmdhole/hole_tcl/reference_bin/` (local
only, never committed; tests that need them skip when absent). Do not use a
previously-installed `~/hole2/exe/hole` as a fidelity reference — an
installed binary can be a stale or locally-patched build — and never
overwrite a user's installed `~/hole2/exe/*`.
