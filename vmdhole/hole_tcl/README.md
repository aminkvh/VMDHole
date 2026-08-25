# hole.tcl — HOLE 2 in pure Tcl (standalone)

A from-source replication of HOLE 2, with no compiled dependency. Standalone so
it can be validated on its own; written with VMDHole in mind so a working
engine can be merged in later as a fallback beside the MOLE one (same `.sph`
record, same profile shape).

Run it:

```sh
tclsh hole.tcl -pdb 1GRM.pdb -rad simple.rad -cvect "0 0 1" \
               -sample 0.5 -endrad 8.0 -csv profile.csv -sph out.sph
```

Defaults match the real binary's own banner: `MCSTEP 1000`, `MCLEN 0.1`,
`kT 0.1`, `sample 0.25`, `endrad 22.0`.

## What works today

Ported by reading the Apache-licensed Fortran directly (`hole2/src`):

- **The RNG, exactly.** This is the piece that makes validation possible at
  all — see below.
- PDB input, `.rad` radius files (the real `VDWR` card — an earlier version of
  this gated on a `RADIUS` card that no real `.rad` file contains, so every
  atom silently used the element fallback; fixed). **CORRECTED, was wrong:**
  `hole::radius_for` used to fall back to `hole::element_radius` (1.85,
  carbon) for ANY atom no `VDWR` rule matched, even when a real radius file
  WAS supplied — a real divergence from the binary, found while picking a
  second validation structure: `tsatr.f` has no such fallback — an atom with
  no matching `VDWR` rule makes the real binary ABORT ("Cannot find vdW
  radius for atom: ..."), not substitute a default (confirmed on a
  K+-containing PDB with `simple.rad`). Fixed: `radius_for` now errors the
  same way when a radius file was supplied and does not cover the atom;
  `element_radius` is only reached when no radius file was supplied at all
  (a separate, still-unaudited case — see that proc's own comment). See
  `vmdhole/hole_tcl/tests/fixtures/1BL8_protein_only.pdb`'s REMARK and
  `vmdhole/hole_tcl/tests/sph_addend_test.sh`'s own abort check for the worked
  example.
- `HOLEEN` (`holeen.f`) — the objective function, energy = −(distance to the
  nearest atom surface), with the three-nearest bookkeeping callers need.
- `HONEWP` (`honewp.f`) — the in-plane perturbation.
- The `HOLCAL` (`holcal.f`) annealing loop, Metropolis acceptance, and the
  linear cooling schedule (`MCKT -= MCKTIN/(0.9·MCSTEP)`, floored at 0).
- `CGUESS`'s `CPOINT` search (`cguess.f`) — CA (or all-atom) centroid, then a
  5-cycle 3×3×3/1 Å hill-climb — used whenever `-cpoint` is not given, in
  place of a plain centroid.
- Axial sampling in both directions with the `ENDRAD` stop, each slice seeded
  from the previous slice's answer; the slice that reaches `ENDRAD` is not
  stored, matching `holcal.f`'s own storage guard.
- `.sph` output, **byte-identical** to `the reference-build hole (see vmdhole/hole_tcl/README.md)`'s own `sphpdb`
  output — not just the stock per-slice `ATOM` records `sph_process` reads,
  but the FULL record: real `STRNOP`/`-STRNON` `IREC` numbering and
  discovery order (NOT sorted by axial coordinate — a plain `lsort` cannot
  reproduce it), the `ADDEND` terminal-sphere grid and its `LAST-REC-END`
  markers (`addend.f`, `hole::addend`), and the duplicate record-0 rewrite at
  the +ve/-ve direction transition (`holcal.f` ~758-806). `hole::write_sph`
  (spherical), `hole::write_capsule_sph` (`QC1`/`QC2` two-atom-per-slice,
  `wpdbsp.f:73-83`) and `hole::write_connolly_sph` (`concal.f`'s own
  multi-record-per-slice write, `concal.f:535-556`) all reproduce it —
  verified `diff -q` clean against `the reference-build hole (see vmdhole/hole_tcl/README.md)` on 1GRM across
  multiple seeds and both an axis-aligned and a non-axis `CVECT`; see
  `vmdhole/hole_tcl/tests/sph_addend_test.sh` for the durable regression and
  `hole::write_sph`'s own header in `hole.tcl` for the exact reproducing
  command. Before this, `write_sph` wrote per-slice records sorted by axial
  coordinate and numbered 0..N-1 positionally — 82 lines vs the real
  binary's 421 on the fixture above — and capsule/CONNOLLY had no `.sph`
  writer at all.
- The 3D surface pipeline: `sph_process.tcl` (`hole::sph_process`, the
  `-sos -dotden D` dot-generation path of `sph_process.f`) and
  `sos_triangle.tcl` (`hole::sos_triangle`, the advancing-front triangulator
  in `sos_triangle.c` — **not** `trisphere.f`, see below). Verified
  **byte-identical, end to end**, against the real binaries on 1GRM at
  dotden 10: `.sph` → `.sos` → `vmd_plot`, `diff -q` clean both ways
  (2657/2657 dots, 2192/2192 after dedup, 4625/4625 raw triangles, 3280/3280
  final triangles). Two single-precision truncations most ports would miss
  had to be replicated exactly — `sphqpu.f`'s dot position (`RVEC3` is
  `REAL`) and `ptgen.f`'s ring-angle ratio (evaluated in single precision
  because `PI` isn't yet part of the expression when it's computed) — see
  `sph_process.tcl`'s header.
- **CONNOLLY** (`concal.f` + `coarea.f`, `vmdhole/hole_tcl/connolly.tcl`) — the flood
  fill that grows a set of probe-radius spheres out from each slice's HOLE
  centre, and the adaptive black/grey/white area estimator that turns them
  into an equivalent radius (Requiv). Verified against the real binary at the
  point level; see `connolly.tcl`'s own header for the full derivation and
  fidelity notes (including the `DIFF3+DIFF3` typo `concal.f` itself has,
  preserved on purpose). **STALE UNTIL THIS REVISION**: this file used to
  list CONNOLLY under "What does not" below as unported/error-only — that was
  wrong by the time `.sph` output (above) started depending on it via
  `hole::write_connolly_sph`; not otherwise re-audited in this revision.
- **CAPSULE** (`vmdhole/hole_tcl/capsule.tcl`) — see that file's own header and the
  measurements later in this document; also not previously listed here.

## What does not

**Turns out not to exist**, for the plain spherical calculation this file
implements — investigated, not ported:

- a post-annealing refinement step. `hsbxmi.f` (401 lines) exists, but it is
  called from exactly one place (`holcal.f:493`) gated `IF (LSPHBX)` — the
  spherebox option, off unless a `SPHPO` card is given. Confirmed against the
  real binary too: it never prints "Applying sd min..." on a plain input. See
  `vmdhole/hole_tcl/refine.tcl` for the full derivation; `hole::refine_slice` exists
  only to error loudly if something ever tries to call it on this path.

So this computes a 2D spherical (and, see `capsule.tcl`, capsule) pore
profile. No number from it should be quoted as a HOLE number without checking
the measurement below still holds for the structure in question.

## The RNG — why this is the interesting part

`holcal.f` calls `DRAND`, which is **not in the Apache source release**: it
comes from a library HOLE links but does not ship. It is recoverable from the
binary — `nm ~/hole2/exe/hole` puts `drand_` at `0x22180`, and disassembling it
shows it seeds once from the `cseed` common block (forcing an odd value) and
then calls `_gfortran_rand`. That is gfortran's `RAND` intrinsic, which
libgfortran implements as the Park–Miller minimal-standard Lehmer generator
with Schrage's method:

```
seed ← 16807 · seed  (mod 2147483647)
value = seed / 2147483647
```

reproduced in `hole::rng`. **A Tcl run can therefore be made to draw the
identical stream as the Fortran one** — which is the only way to distinguish a
real porting error from HOLE's own Monte Carlo noise (measured at 0.0053 Å
between identical reruns of the real binary). Without this the port would be
unvalidatable in principle.

`duran3_` (`0x253f0`), the random unit vector, is verified too — see the
HOLCAL trace section below, which is also where the RNG *sequencing* bug
(3 missing draws once per run, not in this generator at all) was found.

**A silent, unrelated single-precision trap when validating against card
input:** `hole2/src/freda.f`, the card-value parser, has no `IMPLICIT NONE`
and no type declarations, so by Fortran's default implicit typing `FN` (starts
with F) and its digit accumulator `A` are `REAL` — single precision — even
though every physical quantity in the rest of HOLE is `DOUBLE PRECISION`. So
**any decimal number typed into a HOLE control card is silently rounded to
float32 before the run starts** — confirmed by feeding `cpoint -0.0142 -0.0093
4.3477` to the real binary and reading back `DOUBLE(FLOAT32(-0.0142))` etc.
bit-for-bit from an instrumented trace. `hole.f`'s own hardcoded defaults have
the same trap for a different reason: `MCLEN = 0.1` and `MCKTIN = 0.1` are
bare `REAL` literals (no `D0` suffix), so the real binary's actual no-card
default is `DOUBLE(FLOAT32(0.1)) = 0.10000000149011612`, not exact double
`0.1` — `hole::main`'s own `-mclen`/`-mckt` defaults now match this. Neither
of these is a bug in HOLE, and neither is fixable from the Tcl side for a
user-supplied card value (there is no way to bypass FREDA on the real binary)
— but any *test* that compares this port against a card-fed real run and does
not account for it is comparing against a silently-corrupted input, not
against HOLE's real behaviour on the intended numbers.

## Measured against the real binary

1GRM (gramicidin, 1360 atoms), `cvect 0 0 1`, `sample 0.5`, `endrad 8.0`,
`simple.rad`. `~/hole2/exe/hole` is a **stale build** — it is missing a banner
the current source prints, gives a different slice count from a fresh build,
and segfaults on some inputs — so every number below is against a binary
rebuilt from this repo's own `native/stock_build/hole2/src` (`make
clean && make all`; the tree is git-clean before and after). Unpinned
comparisons (`tests/profile_vs_reference.sh`, against whatever seed the real
binary picks for itself) swing between 0.21 and 0.75 Å RMS run to run — that
script is a smoke test, not a metric. Every figure below is a paired-slice
comparison with **both sides pinned to the same seed** (`raseed 1` / `-seed
1`) and an identical explicit `cpoint 0 0 4` (chosen because it, `cvect 0 0
1`, `sample 0.5` and `endrad 8.0` are all exactly representable in the FREDA
single-precision card path above — no decimal-rounding ambiguity to control
for).

| | bottleneck | slices | paired RMS | max |
|---|---|---|---|---|
| `hole.tcl` — spherical | 0.2165 Å @ 6.500 | 40 vs 40 | **1.6e-5 Å** | 1.0e-4 Å |
| `hole.tcl` — capsule, eff.rad | −0.1700 Å @ 2.000 | 31 vs 31 | **0.0003 Å** | 0.0005 Å |
| `hole.tcl` — capsule, cap.rad | | | **0.0002 Å** | 0.0005 Å |

Every slice paired, no isolated outliers, no CLI flags beyond `-cpoint`,
`-cvect`, `-sample`, `-endrad`, `-seed` (defaults for everything else). The
remaining ~0.0005 Å (capsule) is consistent with the real binary's own
**3-decimal-place text output** (`F8.3`/`F12.3` format specifiers) rather than
a residual algorithmic difference — confirmed separately, at full `E24.16`
precision, in the HOLCAL trace below. Reproduced independently on
`raseed 7`/`-seed 7` (different bottleneck location, same agreement).

**Spherical's remaining ~1e-4 Å is text-output rounding too, not algorithm** -
comparing `hole::profile`'s own full-precision return value (not its 4-decimal
`-csv`) against an `E24.16` HOLCAL trace for this exact run gives max
`|Δradius| = 4.4e-16`, `|Δx| = 1.8e-15`, `|Δy| = 4.4e-16`, `|Δz| = 0.0` across
all 40 slices - a few ULPs, i.e. bit-exact. Closing this last gap (from the
0.0152 Å / 0.0950 Å above) needed a **third root cause**, found and fixed by
a step-by-step trace of the plain spherical path (`hole.tcl`'s own `holcal.f`/
`honewp.f`/`holeen.f` triad, done independently of the capsule trace above
and cross-checked against it): `hole::profile`'s reverse-pass restart used
`CPOINT` itself, but `holcal.f`'s own reset (~824-827, "back to original
point") is `LOWCEN(i) = STRCEN(i,0) + SAMPLE*CVECT(i)` with `SAMPLE` already
negated - `STRCEN(,0)` is the **stored, annealed** slice-0 centre, not
`CPOINT`. On 1GRM the anneal moves slice 0's in-plane (x,y) position ~0.2 Å
off `CPOINT` even though `CPOINT` was its own seed, so every negative-direction
slice used to start from a measurably different point than the real search -
explaining why the residual above was one-sided (the positive direction,
which never used this reset, already matched near machine precision; the
negative direction carried the whole 0.0152 Å RMS / 0.0950 Å max). Fixed by
capturing the first stored row's centre into a `slice0` variable during the
forward pass and seeding the reverse pass's `k==0` repositioning from that,
matching `STRCEN(,0)` instead of `CPOINT`. This is scoped to runs with no
`-sph` output: the real binary's `IF (FPDBSP.NE.'NONE')` block (holcal.f
~800-812) resets `LOWCEN` a *second* way before this one, from
`STRCEN(,STRNOP-1)` ("the one before as the last is sometimes wonky"), which
`hole_tcl` does not implement - not exercised here since no `-sph` was given
on either side, but real for anyone who does pass one.

An earlier revision of this file reported capsule figures of eff.rad max
2.68 Å / RMS 0.58 Å, cap.rad max 0.71 Å / RMS 0.23 Å, and a 39-vs-36 slice
mismatch, with spherical showing three isolated points diverging by up to
3.09 Å despite 38 of 41 paired points agreeing to 0.0066 Å. **Those numbers
were real** (reproducible on that revision of the code) but were misdiagnosed
as search-algorithm or floating-point noise; the two real causes, found by a
step-by-step trace (see below), are fixed as of this revision.

### The HOLCAL trace: finding the real divergence

The method: add `WRITE(NOUT,...)` calls at full `E24.16` precision to a
**local, uncommitted** copy of `holcal.f`/`honewp.f` (never `LDBUG`, which is
hardcoded `.FALSE.` with no card to set it — these were plain unconditional
writes), rebuild, run, diff against an identically-instrumented Tcl trace
harness that calls the same `hole::hcapen`/`hole::honewp`/`hole::rng::*`
primitives capsule.tcl uses, then `git checkout` the Fortran source back to
pristine and rebuild clean again. A previous round of trying to find this bug
by reasoning about the code and by comparing final profiles could not localise
it; comparing every Metropolis step let the trace name the exact first
divergent value instead of guessing.

**Root cause 1 (the dominant one): a missing RNG draw.** `hole.f:548-549`,
called once per run, unconditionally, before `HOLCAL` starts and *after*
`CVECT` is unitised:

```fortran
C New feature 24 October 1993 write out seed integer used
C find one random 3d vector to kick things off
      CALL dURAN3( DISP)
```

`DISP` is never used for anything — the call exists purely to consume 3
`DRAND` draws. This port never made it, so every real HOLE run's RNG stream
was offset by exactly 3 draws relative to this port's, from the very first
`HONEWP` call onward. Found by: `hole::rng::seed_like_drand(1)` followed by 8
`hole::rng::rand()` calls reproduces an isolated ground-truth probe (`FSEED=1`,
raw `DRAND` calls, linked against the real `hole.a`) bit-for-bit — so the RNG
*model* was never wrong — but the real run's *first* `HONEWP` call's unit
vector matched `duran3` applied to probe draws 4–6, not 1–3, meaning 3 draws
were consumed by something before `HONEWP` got a turn. That something is
`hole.f:549`. Fixed as `hole::rng::kick_off`, called once from `hole::main`
right after `hole::rng::seed_like_drand`.

**Root cause 2: `MCLEN`/`MCKTIN`'s real no-card default is not exact double
0.1.** See the single-precision-card-trap note above — `hole.f`'s `MCLEN =
0.1` / `MCKTIN = 0.1` are bare `REAL` literals, so the real binary's default is
`DOUBLE(FLOAT32(0.1)) = 0.10000000149011612`. Found the same way: with root
cause 1 fixed, `HONEWP`'s step-2 `NEWCEN` still disagreed by ~1e-8 relative
even with the RNG stream and `CPOINT` exactly aligned; the trace showed
`MCLEN` itself printing as `0.1000000014901161` on the Fortran side, not
`0.1000000000000000`. `hole::main`'s `-mclen`/`-mckt` defaults now match.

Two smaller, source-confirmed (not trace-required) fixes landed alongside
these: `HONEWP`'s normalise-then-scale is now done as three separate
per-component divisions followed by two separate multiplications, matching
`honewp.f`'s actual `dUVEC2` call and `DISP(i)=DISP(i)*MCLEN*WORK` — not the
single combined scalar this port used before, which rounds differently in the
last bit — and the Metropolis `PROB = EXP(-HGOOD/MCKT)` no longer has an
`arg < -700 ? 0.0 : exp(arg)` clamp that doesn't exist in the source (Tcl's
`exp()` underflows to `0.0` on its own, same as Fortran's, confirmed directly:
`exp(-720.0)` → `2.03e-313`, no domain error).

With both root causes fixed, a full step-by-step trace of a 42-slice, 1000-MC-
step-per-slice capsule search (30-step reduced version for tractability) shows
**zero** tag/branch mismatches — every single accept/reject, collapse-vs-keep-
moved, and storage decision matches the real binary exactly — with only the
~1e-9-relative residual described above (traced to `HCAPEN`'s own per-atom
loop; ruled out atom-table precision by comparing every coordinate/radius the
real `TSATR`/`TSRADR` readers produce against `hole::read_pdb`'s, bit-for-bit
identical). This is a materially different, much smaller kind of residual than
the isolated multi-Å outliers reported previously, and it never once flipped a
branch decision across the whole traced run.

The **39-vs-36 (and 41-vs-38) slice-count mismatches investigated previously
were not a separate bug** — the storage guard (`IF (-LOWENG.LT.ENDRAD)`,
matching `capsule.tcl`'s own `if {$effrad >= $o(-endrad)} { break }`) was
already correct; the count differences were a symptom of the RNG/MCLEN
misalignment producing a genuinely different search trajectory. With the two
fixes above, slice counts now match exactly (checked on two independent
seeds).

**Root cause 3 (closing the remaining ~1e-9-relative residual): `CAPSULE_PI`
was true double pi, not the real binary's `DOUBLE(FLOAT32(2.*ASIN(1.)))`.**
`hole.f:363` computes `PI = 2.*ASIN(1.)` — bare `REAL` literals, so `ASIN`'s
argument/result and the multiply are done in float32 and only widened to
double on assignment, the same silent-single-precision trap as `MCLEN`
above. This `PI` is passed into `HCAPEN` as an argument and used **only** in
the area transform (`ENERGY = PI*ENERGY**2 + 2*ENERGY*DCENT`, then
`SQRT(.../PI)`) — `CAPRAD`, the raw pre-transform radius, never touches it.
Confirmed against the *actual instrumented binary* (not a proxy): adding a
`TRANSFER`-to-`INTEGER*8` `WRITE` right after `hole.f:363` in a scratch copy
of the source tree and rebuilding prints hex `400921FB60000000` for `PI` at
runtime, bit-identical to Tcl's `binary scan [binary format f [expr
{2.0*asin(1.0)}]] f y`, vs. `400921FB54442D18` for true double pi.
`hcapen.f:197`'s `PI*ENERGY**2` also associates as `PI*(ENERGY*ENERGY)`, not
`(PI*ENERGY)*ENERGY` — Fortran `**` binds tighter than `*` — a second,
smaller association-order fix applied at both `ENERGY`/`DAT2` sites.

Measured with a **full per-slice E24.16 trace of the whole two-centre
search's stored state** (`holcal.f`'s `STRNOP`/`STRNON` storage sites
instrumented to dump `LOWCEN`, `LOWLVC`, `-LOWENG`, `LOWBRD`; diffed against
a Tcl harness replaying `hole::anneal_slice_capsule` in discovery order —
not the 4-decimal CSV or the 3-decimal printed table used for the table
above). 1GRM, `cpoint 0 0 4`, `cvect 0 0 1`, `sample 0.5`, `endrad 8.0`:

| seed | slices | without PI fix: max|Δeff.rad| | with PI fix: max|Δeff.rad| | max|Δcentre/LVC/dcent| (either way) |
|---|---|---|---|---|
| `raseed 1` | 31 = 31 | 9.405e-08 Å | 8.882e-16 Å | ≤3.6e-15 Å |
| `raseed 7` | 36 = 36 | (not separately measured) | 4.441e-16 Å | ≤3.6e-15 Å |

Index sets matched exactly on both seeds (e.g. seed 1: `-18..12`, no gaps).
**The two-centre separation (`dcent`) and both centres were already
bit-exact — at the 1e-15 floor — even *before* the PI fix**: only the
PI-dependent area transform was off, by ~1.4e-8 relative, five orders of
magnitude below the 0.0005 Å print quantum in the table above, which is why
it was invisible at CLI/CSV precision and is the direct proof that the
0.0003/0.0005 Å capsule figures are formatting, not algorithm. The
association-order and `dcent`-threshold (`1E-09` → `DOUBLE(FLOAT32(1e-9))`)
edits are correct per source but **measurably inert** on both traced seeds —
identical output whether applied or not — and the `dcent` one is provably so
regardless of fixture: the two candidate thresholds differ by ~2.8e-17, and
every `dcent==0.0` case this port's own search produces (the
collapse-to-midpoint reseed calling `hcapen` with centre exactly equal to
seccen) is already far below either threshold identically. Kept for source
fidelity, not because either changed a measured value here.

This also **settles the "two-centre separation diverges" question directly**
rather than by elimination: that language is from commits `73cb4ab7` and
`63dcea33`, both measured *before* the RNG-kick-off/MCLEN fixes in
`02d26f00` — i.e. on a desynchronized RNG stream. At the current revision
the search state (not just the reported radius) is bit-exact. See
`capsule.tcl`'s `CAPSULE_PI` definition for the full reproducing recipe.

## Generalisation beyond 1GRM/z (capsule)

**Done**, not speculative — `vmdhole/hole_tcl/tests/capsule_generalization_test.sh`
(print-precision regression guard, runs as part of `run_all.sh`) and, for the
bit-exact claim, the same full E24.16 trace method as "Root cause 3" above,
run on two more configurations:

- **Non-axis `CVECT`**: 1GRM, `cvect 1 1 1` (`cpoint 0 0 4`, `sample 0.5`,
  `endrad 8.0`, `raseed 1`) — 33 = 33 slices, max|Δeff.rad| 4.4e-16 Å,
  max|Δcentre/LVC/dcent| 3.6e-15 Å. Bit-exact, same floor as the axis-aligned
  case.
- **A second, independent structure**: KcsA (`1BL8`, 4 chains, 2820 atoms —
  ~2x 1GRM's atom count), protein atoms only (`vmdhole/hole_tcl/tests/fixtures/
  1BL8_protein_only.pdb` — the original's K+/HOH HETATM records make the
  *real* binary abort with "Cannot find vdW radius for atom: K", since
  `tsatr.f` has no element fallback of its own; see the fixture's REMARK
  header), `cvect 0 0 1`, **no `-cpoint`** (CGUESS on both sides — this also
  exercises CGUESS on a structure it was never checked against before),
  `sample 0.5`, `endrad 10.0`, `raseed 1` — 115 = 115 slices, max|Δeff.rad|
  8.9e-16 Å, max|Δcentre/LVC/dcent| 4.3e-14 Å (absolute; the coordinates
  themselves are ~10-100 Å here, so still ~1e-16 relative). CGUESS produced
  the identical `CPOINT` on both sides (`74.2201881443... 26.5988969072...
  24.9731855670...`, matching the real binary's printed `74.2202 26.5989
  24.9732` to the digits it shows).

Both traces used the same instrumented-tree recipe as "Root cause 3", with
a harness generalized to accept `cpoint`/`cvect`/`sample`/`endrad` as
arguments (`hole::cguess_cpoint` when no `cpoint` is given) rather than the
hardcoded 1GRM values — not checked in (scratch-only), reproducible with the
same recipe.

So the capsule bit-exactness claim is **not scoped to axis-aligned CVECT or
to 1GRM** — both a genuinely different structure and a genuinely
non-axis-aligned direction land at the same ~1e-15/1e-16 floor, and 1BL8's
own bottleneck (eff.rad down to -0.095 Å with mid-points, i.e. the fitted
capsule genuinely overlaps atom vdW spheres there) is a materially tighter
constriction than 1GRM's. What was *not* tried: a slice where the two-centre
separation (`DCENT`) itself gets genuinely small - the smallest traced value
across both new configurations was 0.034 Å (1BL8) - three-plus orders of
magnitude above the `1E-09` threshold in `hcapen.f:101`; every `DCENT==0.0`
this port's own search produces comes from the exact degenerate
collapse-to-midpoint reseed (see "Root cause 3" above), not from the
annealing search organically converging to near-zero separation. That region
of `HCAPEN`'s branch remains unexercised by any measurement in this file.

## Next steps, in order

1. `sph_process.tcl`/`sos_triangle.tcl` are ~150-200x slower than the real
   binaries (same O(n) complexity, Tcl interpreter overhead) — measured
   ~8.8 s vs ~54 ms end to end on 1GRM's 199-sphere `.sph` at dotden 10; a
   trajectory surface is not yet practical without the binaries. CONNOLLY's
   own `.sph` (many records per slice, not one) makes this worse still —
   feeding a CONNOLLY `.sph` through `hole::sph_process` was not timed to
   completion in the session that added `hole::write_connolly_sph` (still
   running past 2 minutes on a small fixture where the real binary takes
   under a second); use the real `sph_process`/`sos_triangle` binaries for a
   CONNOLLY-mode `.sph`, not this port's own triangulator, until this is
   revisited.
2. `hole::sph_process` (the `.sph` READER, not this file's own writers) does
   not understand the capsule `QC1`/`QC2` record — `hole::write_capsule_sph`
   produces a `.sph` `the reference build's sph_process` renders correctly, but
   `hole::sph_process` itself raises "CAPSULE .sph records (QC1/QC2) are not
   ported" on the same file. Not fixed in the session that found it (a
   different subsystem to the write side this session was scoped to) — see
   `hole::_sph_read` in `sph_process.tcl`.

## Licensing

HOLE 2 is Apache License 2.0 (`hole2/LICENSE`), so both this derivation and
redistribution of the original binaries are permitted. Worth remembering that
the second option solves the "user has no working HOLE" problem for all three
binaries at a fraction of the cost — see
`vmdhole/NOTES/pure-tcl-hole-feasibility.md`.


---

## Reference binaries (formerly `the reference build's `)

`hole`, `sph_process`, `sos_triangle` built from
`native/stock_build/hole2/src` (upstream `osmart/hole2` commit
`a8eaf6121ba66625446933f4acd7d6aa336dbb47`), pristine - no local patches -
using the documented recipe (`hole2/INSTALL.md`):

```sh
cd native/stock_build/hole2/src
source ../source.apache      # sets HoleVersion=2.3.1 etc. so the banner isn't the default
make ../exe/hole ../exe/sph_process ../exe/sos_triangle
```

Built in an isolated private copy of the source tree (not in place in
`native/stock_build/`, which is shared with other concurrent work and
was observed mid-build - a stray `vertim.o`/`hole.a` race is easy to hit
there), then copied here. gfortran 13.3.0 (Ubuntu 24.04), gcc default `cc`.

## `~/hole2/exe/hole` must not be used as a fidelity reference

Confirmed concretely (1GRM, `simple.rad`, `cvect 0 0 1`, `sample 0.5`,
`endrad 8.0`):

1. **Banner differs**: `~/hole2/exe/hole` prints `HOLE release SOURCE
   DISTRIBUTION (?)` / `3rd party build...`; this binary prints `HOLE release
   2.3.1 (01 February 2023)` / `usage subject to Apache License, Version
   2.0` - expected, since `~/hole2/exe/hole` was not built via
   `source ../source.apache`.
2. **`CONNOLLY` reliably SEGFAULTs** on `~/hole2/exe/hole` - reproduced at
   `raseed 1, 2, 3, 42` and with no `raseed` card at all, every time. This
   binary never crashes on the same inputs.
3. **Root cause identified, not just observed**: `~/hole2/exe/hole` is
   dynamically linked against `libgomp.so.1` and defines
   `holcal_._omp_fn.0` (`nm`/`ldd`) - i.e. it is a build of
   `native/connolly_patches/` (`holcal_par.f` compiled `-fopenmp`,
   swapped in for stock `holcal.f`), the project's own experimental
   parallel-CONNOLLY acceleration patch, not stock/Apache HOLE. Applying
   `connolly_patches/apply_patches.py` to *this* pristine tree and rebuilding
   reproduces the identical symbol signature (`libgomp`, `holcal_._omp_fn.0`,
   `GOMP_parallel` etc.) **and** the identical CONNOLLY+RASEED SIGSEGV. So
   the crash is very likely a threading bug in that experimental patch set,
   not evidence against the *stock* HOLE algorithm or against this port.
4. **Slice-count claim does not survive seed control.** At a pinned seed
   (`raseed 1`) `~/hole2/exe/hole` and this build produce **byte-identical**
   stdout (42/42 slices, every column) apart from the 4 banner lines.
   Unseeded, both binaries independently vary 39-44 slices run to run with
   >1s between runs (ordinary HOLE Monte Carlo noise, not a build artefact) -
   the ranges fully overlap, so a single unseeded run of each binary is not
   evidence of a real difference. **Retracted**, with numbers, as a symptom
   of staleness.
5. **Spherical and CAPSULE paths are numerically unaffected** by whichever
   patch set produced `~/hole2/exe/hole`: at `raseed 1`, both plain-spherical
   and `CAPSUL` runs are byte-identical to this build apart from the banner.
   Only `CONNOLLY` diverges (by crashing). This means the previously-recorded
   spherical (38/41 @ 0.0066 A) and capsule (RMS 0.0004 A) fidelity numbers in
   `vmdhole/hole_tcl/README.md`, even though measured against `~/hole2/exe/hole`,
   were not measuring an artefact of that binary's patches - see
   `vmdhole/hole_tcl/tests/` and the delegate report for the re-verification.

Do not overwrite `~/hole2/exe/*` (the user's install) or treat it as ground
truth for CONNOLLY. Use this directory instead.


(The binaries themselves are built locally per the steps above and are
never committed; tests that need them skip when they are absent.)
