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

## What is ported

Everything the plain spherical, CAPSULE and CONNOLLY calculations need, ported
by reading the Apache-licensed Fortran directly: PDB and `.rad` input (an atom
no `VDWR` rule covers aborts the run, as in the real binary), the annealing
loop with Metropolis acceptance and cooling schedule, the `CPOINT` guess,
bidirectional axial sampling with the `ENDRAD` stop, and full `.sph` output —
**byte-identical** to the reference binaries' own `sphpdb` records, including
`IREC` numbering/discovery order and the `ADDEND` terminal grid
(`tests/sph_addend_test.sh`). The 3D surface pipeline (`sph_process.tcl`,
`sos_triangle.tcl`) is byte-identical end to end on the tested fixtures. A
post-annealing refinement step is deliberately absent: it only exists upstream
under the spherebox option (see `refine.tcl`).

Per-file derivation and fidelity notes live in each source file's own header
(`hole.tcl`, `capsule.tcl`, `connolly.tcl`, `sph_process.tcl`).

## Validation (the HOLCAL trace)

The port draws the **identical random stream** as the real binary:
`DRAND` is gfortran's `RAND` intrinsic (Park–Miller/Schrage), reproduced in
`hole::rng`. That is what makes exact validation possible — HOLE's own
Monte Carlo noise between identical reruns is larger than any porting error
of interest.

Bit-exactness was established by step-by-step full-precision (`E24.16`)
tracing against an instrumented local build, which found and fixed three
divergences: a missing once-per-run RNG kick-off draw (root cause 1), the
real binary's float32 no-card defaults (root cause 2), and `CAPSULE_PI` being
`DOUBLE(FLOAT32(2.*ASIN(1.)))` rather than true double pi (root cause 3 — see
`capsule.tcl`'s `CAPSULE_PI` header). At pinned seed and `cpoint`, spherical
and capsule agree with the real binary to ~1e-15 Å across every slice, on two
structures (1GRM, 1BL8) and a non-axis `CVECT`; the ≤0.0005 Å residuals in
printed output are the binary's 3-decimal text format, not algorithm.
`tests/capsule_generalization_test.sh` is the regression guard.

**When validating against card input, note** that HOLE's card parser
(`freda.f`) silently rounds every decimal card value to float32 before the run
starts. Choose exactly-representable values (or match the rounding) or the
comparison is against corrupted input. Unpinned comparisons
(`tests/profile_vs_reference.sh`) swing ~0.2–0.8 Å RMS run to run from HOLE's
own Monte Carlo noise — that script is a smoke test, not a metric.

## Next steps

1. `sph_process.tcl`/`sos_triangle.tcl` are ~150–200× slower than the real
   binaries (Tcl interpreter overhead) — use the real binaries for trajectory
   surfaces and for CONNOLLY-mode `.sph` files.
2. `hole::sph_process` (the `.sph` **reader**) does not understand the capsule
   `QC1`/`QC2` record; the real `sph_process` renders those files correctly.
   See `hole::_sph_read` in `sph_process.tcl`.

## Licensing

HOLE 2 is Apache License 2.0 (`hole2/LICENSE`), so both this derivation and
redistribution of the original binaries are permitted.

## Reference binaries

Fidelity comparisons use `hole`, `sph_process`, `sos_triangle` built from
`native/stock_build/hole2/src` (upstream `osmart/hole2` commit
`a8eaf6121ba66625446933f4acd7d6aa336dbb47`), pristine — no local patches:

```sh
cd native/stock_build/hole2/src
source ../source.apache      # sets HoleVersion=2.3.1 etc.
make ../exe/hole ../exe/sph_process ../exe/sos_triangle
```

then copy the three binaries into `vmdhole/hole_tcl/reference_bin/` (local
only, never committed; tests that need them skip when absent).

**Do not use a previously-installed `~/hole2/exe/hole` as a fidelity
reference** — an installed binary can be a stale or locally-patched build.
Make fidelity claims against a pristine build, and never overwrite a user's
installed `~/hole2/exe/*`.
