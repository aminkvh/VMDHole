# hole.tcl — HOLE 2 in pure Tcl (standalone)

A from-source replication of HOLE 2 with no compiled dependency. Standalone so
it can be validated on its own; the same engine is inlined into `vmdhole.tcl`
as the plugin's pure-Tcl fallback (same `.sph` record, same profile shape).

Run it:

```sh
tclsh hole.tcl -pdb 1GRM.pdb -rad simple.rad -cvect "0 0 1" \
               -sample 0.5 -endrad 8.0 -csv profile.csv -sph out.sph
```

Defaults match the real binary's own banner: `MCSTEP 1000`, `MCLEN 0.1`,
`kT 0.1`, `sample 0.25`, `endrad 22.0`.

## What is ported

Ported by reading the Apache-licensed Fortran directly (`hole2/src`):

- **The RNG, exactly** — the piece that makes validation possible at all; see
  below.
- PDB input and `.rad` radius files (the real `VDWR` card). As in the real
  binary (`tsatr.f`), an atom no `VDWR` rule covers makes the run **abort**,
  not fall back to a default radius; the element fallback is only reached when
  no radius file was supplied at all.
- `HOLEEN` (`holeen.f`) — the objective function, with the three-nearest
  bookkeeping callers need.
- `HONEWP` (`honewp.f`) — the in-plane perturbation.
- The `HOLCAL` (`holcal.f`) annealing loop, Metropolis acceptance, and the
  linear cooling schedule.
- `CGUESS`'s `CPOINT` search (`cguess.f`), used whenever `-cpoint` is not
  given.
- Axial sampling in both directions with the `ENDRAD` stop, matching
  `holcal.f`'s own storage guard.
- `.sph` output, **byte-identical** to the reference binaries' own `sphpdb`
  output — the full record: real `STRNOP`/`-STRNON` `IREC` numbering and
  discovery order, the `ADDEND` terminal-sphere grid with its `LAST-REC-END`
  markers, and the duplicate record-0 rewrite at the direction transition.
  `hole::write_sph` (spherical), `hole::write_capsule_sph` (`QC1`/`QC2`) and
  `hole::write_connolly_sph` (multi-record-per-slice, `concal.f:535-556`) all
  reproduce it; `tests/sph_addend_test.sh` is the durable regression.
- The 3D surface pipeline: `sph_process.tcl` and `sos_triangle.tcl`, verified
  **byte-identical end to end** against the real binaries on 1GRM at dotden 10
  (`.sph` → `.sos` → `vmd_plot`, both ways). Two single-precision truncations
  had to be replicated exactly — see `sph_process.tcl`'s header.
- **CONNOLLY** (`connolly.tcl`) — the probe-sphere flood fill and the adaptive
  area estimator (Requiv); see that file's header for the derivation and
  fidelity notes.
- **CAPSULE** (`capsule.tcl`) — see that file's header and the measurements
  below.

## What is not

A post-annealing refinement step **turns out not to exist** for the plain
spherical calculation: `hsbxmi.f` is only called under the spherebox option
(`SPHPO` card), confirmed against the real binary. `hole::refine_slice` exists
only to error loudly if something ever calls it on this path (see
`refine.tcl`).

## The RNG — why validation is possible

`holcal.f` calls `DRAND`, which is not in the Apache source release; it is
recoverable from the binary. Disassembly shows it seeds once from the `cseed`
common block (forcing an odd value) and then calls gfortran's `RAND`
intrinsic — the Park–Miller minimal-standard Lehmer generator with Schrage's
method:

```
seed ← 16807 · seed  (mod 2147483647)
value = seed / 2147483647
```

reproduced in `hole::rng`. A Tcl run can therefore draw the **identical
stream** as the Fortran one — the only way to distinguish a real porting
error from HOLE's own Monte Carlo noise (measured at 0.0053 Å between
identical reruns of the real binary).

**A single-precision trap when validating against card input:** `freda.f`,
HOLE's card parser, accumulates values in `REAL` (implicit typing), so **any
decimal number typed into a HOLE control card is silently rounded to float32
before the run starts**. `hole.f`'s own hardcoded defaults have the same trap
(`MCLEN = 0.1` etc. are bare `REAL` literals, so the real default is
`DOUBLE(FLOAT32(0.1))`, not exact double 0.1) — this port's defaults match
the real binary's actual values. Any comparison against a card-fed real run
must account for this or it compares against silently-corrupted input.

## Measured against the real binary

1GRM (gramicidin, 1360 atoms), `cvect 0 0 1`, `sample 0.5`, `endrad 8.0`,
`simple.rad`, both sides pinned to the same seed (`raseed 1` / `-seed 1`) and
an explicit `cpoint 0 0 4` (all values exactly representable in the FREDA
float32 card path, so no rounding ambiguity):

| | bottleneck | slices | paired RMS | max |
|---|---|---|---|---|
| `hole.tcl` — spherical | 0.2165 Å @ 6.500 | 40 vs 40 | **1.6e-5 Å** | 1.0e-4 Å |
| `hole.tcl` — capsule, eff.rad | −0.1700 Å @ 2.000 | 31 vs 31 | **0.0003 Å** | 0.0005 Å |
| `hole.tcl` — capsule, cap.rad | | | **0.0002 Å** | 0.0005 Å |

(Unpinned comparisons — `tests/profile_vs_reference.sh`, against whatever seed
the real binary picks — swing 0.21–0.75 Å RMS run to run from HOLE's own Monte
Carlo noise; that script is a smoke test, not a metric.)

The residuals in that table are the real binary's **3-decimal text output**,
not algorithm: compared at full `E24.16` precision (see the trace below), both
spherical and capsule agree to ~1e-15/1e-16 Å — bit-exact. Reproduced
independently on a second seed.

## The HOLCAL trace (root causes 1–3)

Final profiles could not localise the early multi-Å divergences, so the port
was validated by step-by-step tracing instead: `WRITE` calls at `E24.16`
precision added to a local, uncommitted copy of the Fortran, diffed against an
identically-instrumented Tcl harness, then the source restored. Comparing
every Metropolis step names the exact first divergent value instead of
guessing. Three root causes were found and fixed this way:

1. **A missing RNG draw.** `hole.f:548-549` unconditionally calls
   `dURAN3(DISP)` once per run before `HOLCAL` starts, purely consuming 3
   `DRAND` draws (`DISP` is never used). Without it every real run's stream
   was offset by exactly 3 draws. Fixed as `hole::rng::kick_off`.
2. **The real no-card defaults are not exact double 0.1.** See the
   single-precision trap above; `-mclen`/`-mckt` defaults now match
   `DOUBLE(FLOAT32(0.1))`.
3. **`CAPSULE_PI` is the real binary's `DOUBLE(FLOAT32(2.*ASIN(1.)))`, not
   true double pi.** `hole.f:363` computes `PI = 2.*ASIN(1.)` in float32
   (bare `REAL` literals) — confirmed by instrumenting the binary (hex
   `400921FB60000000` at runtime). This `PI` only enters the capsule area
   transform; with the fix the remaining ~1e-8-relative residual closes to
   ~1e-16. `hcapen.f`'s `PI*ENERGY**2` association order is also matched. See
   `capsule.tcl`'s `CAPSULE_PI` definition for the reproducing recipe.

With all three fixed, a full trace of a 42-slice capsule search shows **zero**
accept/reject or storage decision mismatches against the real binary, and the
reverse-pass restart matches `holcal.f`'s own reset (`STRCEN(,0)`, the stored
annealed slice-0 centre — not `CPOINT`). One caveat: the real binary's
`IF (FPDBSP.NE.'NONE')` block resets `LOWCEN` a second way when a `-sph`
output is requested (`STRCEN(,STRNOP-1)`); that variant is not implemented.

## Generalisation beyond 1GRM/z (capsule)

Verified with the same full `E24.16` trace method as Root cause 3, not
assumed (`tests/capsule_generalization_test.sh` is the print-precision
regression guard):

- **Non-axis `CVECT`**: 1GRM, `cvect 1 1 1` — 33 = 33 slices,
  max|Δeff.rad| 4.4e-16 Å.
- **A second structure**: KcsA (`1BL8`, protein atoms only — the K+/HOH
  HETATM records make the *real* binary abort, see the fixture's REMARK
  header), no `-cpoint` (exercising CGUESS on both sides) — 115 = 115 slices,
  max|Δeff.rad| 8.9e-16 Å, identical CGUESS `CPOINT` on both sides.

So the bit-exactness claim is not scoped to axis-aligned CVECT or to 1GRM.
Not exercised by any measurement here: a slice where the two-centre
separation (`DCENT`) organically becomes near-zero (smallest traced value
0.034 Å, far above the `1E-09` threshold in `hcapen.f:101`).

## Next steps, in order

1. `sph_process.tcl`/`sos_triangle.tcl` are ~150–200× slower than the real
   binaries (same O(n) complexity, Tcl interpreter overhead) — measured
   ~8.8 s vs ~54 ms end to end on 1GRM's 199-sphere `.sph` at dotden 10. A
   trajectory surface is not practical without the binaries, and a
   CONNOLLY-mode `.sph` (many records per slice) is worse still — use the
   real `sph_process`/`sos_triangle` binaries for those.
2. `hole::sph_process` (the `.sph` **reader**) does not understand the
   capsule `QC1`/`QC2` record — `hole::write_capsule_sph` produces a `.sph`
   the real `sph_process` renders correctly, but this port's own reader
   raises "CAPSULE .sph records (QC1/QC2) are not ported" on the same file.
   See `hole::_sph_read` in `sph_process.tcl`.

## Licensing

HOLE 2 is Apache License 2.0 (`hole2/LICENSE`), so both this derivation and
redistribution of the original binaries are permitted.

---

## Reference binaries

Fidelity comparisons use `hole`, `sph_process`, `sos_triangle` built from
`native/stock_build/hole2/src` (upstream `osmart/hole2` commit
`a8eaf6121ba66625446933f4acd7d6aa336dbb47`), pristine — no local patches —
using the documented recipe:

```sh
cd native/stock_build/hole2/src
source ../source.apache      # sets HoleVersion=2.3.1 etc. so the banner isn't the default
make ../exe/hole ../exe/sph_process ../exe/sos_triangle
```

then copy the three binaries into `vmdhole/hole_tcl/reference_bin/` (local
only, never committed; tests that need them skip when they are absent). Build
in a private copy of the source tree — the shared one can be mid-build.

**Do not use a previously-installed `~/hole2/exe/hole` as a fidelity
reference.** A locally-installed binary can be a stale or locally-patched
build: the one found here printed a different banner, was dynamically linked
against `libgomp` (i.e. built from the experimental parallel-CONNOLLY patch
set, not stock), and reliably segfaulted on CONNOLLY runs. At a pinned seed
its spherical and capsule output was still byte-identical to the pristine
build, so numbers measured against it remain valid for those modes — but
fidelity claims should be made against a pristine build, and a user's
installed `~/hole2/exe/*` should never be overwritten.
