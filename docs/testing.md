# Testing

VMDHole's tests live in three tiers, plus a set of native byte-identity
verifications. Each tier exists for a different question, needs different
prerequisites, and is safe to run at a different frequency. Everything here
follows two house rules:

1. **Every test is verified to go RED on the real defect before being
   trusted.** A green test that could not have caught the bug it is named
   after is worse than no test. Unit-test headers record the defect each
   guards (most also record how the RED run was produced); for the main
   suite the per-group defect catalogue and the RED-verification statement
   live in `run_tests.sh`'s own header.
2. **A test that genuinely cannot run prints `SKIP:` at column 0 and exits
   0.** The release gate (`VMDHOLE_RELEASE=1`) anchors on that form to tell
   "checked nothing" from a real pass;
   `tests/unit/test_skip_message_contract.sh` polices it.

## The tiers at a glance

| Tier | Where | Question it answers | Needs | Run it with |
|---|---|---|---|---|
| Unit regression | `tests/unit/` | one defect each, in isolation | `sh`, a C compiler, `tclsh`; Tk + an X display for the GUI smoke test | `sh tests/unit/run_unit_tests.sh` |
| Main suite | `vmdhole/tests/` | the plugin, end to end — 20 groups | VMD for most groups; locally built engines; reference HOLE for parity groups | `vmdhole/tests/run_tests.sh` |
| Pure-Tcl engine | `vmdhole/hole_tcl/tests/` | the Tcl HOLE engine vs a reference binary | a reference `hole` build | `vmdhole/hole_tcl/tests/run_all.sh` (`profile_vs_reference.sh <pdb>` is the by-hand comparison harness) |
| Native verifications | `native/` | byte-identity of the accelerated binaries vs stock HOLE | stock `hole2` tree (source and/or binaries) | `native/verify.sh`, `native/connolly_patches/test_hcapen_cache.sh` |

CI (`.github/workflows/tests.yml`) runs the unit suite under `xvfb-run`
(so the GUI smoke test executes headlessly), plus the C builds, Tcl lint,
release-integrity, and a hand-picked set of VMD-free main-suite groups.

## Environment knobs

| Variable | Effect |
|---|---|
| `VMDHOLE_HOLE_EXE_DIR` | directory searched FIRST for `hole` and its siblings — honoured by the plugin's own discovery (`find_hole_exe`) and by the reference-dependent tests. Point it at `native/build` on a tree that builds there. |
| `VMDHOLE_CONFIG_FILE` | overrides `~/.vmdhole_config`. `run_tests.sh` exports a per-run temp file automatically, so the suite never reads or rewrites the user's real config. |
| `VMD` / `VMD_BIN` | which VMD binary the wrappers drive (`VMD_BIN` for `test_gui_reachable.sh` and `test_adapter_schema.sh`, `VMD` elsewhere). |
| `VMDHOLE_RELEASE=1` | a skipped group becomes a FAILURE — required before tagging. |
| `EXE`, `STOCK`, `SRC` | `test_accel_parity.sh`'s accelerated dir, stock `sph_process`, and patched `hole2/src` tree. |
| `GUI_TEST_*` | `test_gui_reachable.sh`'s fixture overrides (PDB, HET residue, selection, start point, engine). |

## Unit regression tests (`tests/unit/`)

Small and self-contained: no VMD, no trajectory data, no gitignored fixture
corpus, no network. The whole directory is new relative to VMDHole 1.0.0;
tests marked **(review)** were added by the review-and-harden pass that also
produced the fixes they guard.

| Test | Guards |
|---|---|
| `test_sos_nonfinite_input.sh` | signed-integer overflow (UB) in `sos_triangle_fast`'s spatial grids on a `.sos` carrying non-finite coordinates |
| `test_sos_short_record_normals.sh` **(review)** | a short `.sos` record inheriting the previous record's normal (or an indeterminate stack value on record one) |
| `test_skip_message_contract.sh` | group-level skips written in a form the release gate cannot see |
| `test_mole_cli_validation.sh` | the tunnel engine's CLI: negative `--cover` recursion, out-of-range weight, flags parsed as positionals, non-finite/junk numeric argv **(review: the nan/inf cases)** |
| `test_mole_coord_validation.sh` | `nan`/`inf`/junk atom coordinates producing a plausible tunnel set with exit 0 |
| `test_mole_depth_bound.sh` | `remove_shallow` iterating the raw unclamped MinDepth (~12 h at `INT_MAX`) |
| `test_preset_exec_paths.sh` | `load_config` discarding engine paths a batch script set first |
| `test_headless_scratch_reclaim.sh` | `/dev/shm` scratch reclaim reachable only through the GUI |
| `test_headless_run_guards.sh` **(review)** | tunnel start-point validation and shell quoting (nan/inf included), busy/`_end_calc` restoration across a throw, the shared atomselect's lifetime in `run_analysis` |
| `test_tsv_reader_parity.sh` | the threaded TSV reader skipping `_resolve_conn_radii` |
| `test_tsv_publish_on_failure.sh` **(review)** | a failed profile parse truncating the good `hole_profile.tsv` beside it — both writers publish by rename only after a parse that produced rows |
| `test_gui_smoke.sh` **(review)** | the GUI itself without VMD: `vmdhole.tcl` sourced under plain tclsh+Tk with VMD stubbed, the real widget tree built by the real `show_gui`, and scripted user actions asserting the close path, the busy guards, the deleted-molecule dialog class, tunnel gear-popup route pinning, and nan-tolerant option fields. Skips without Tk or a display; CI runs it under xvfb. |

`run_unit_tests.sh` globs `test_*.sh`, so a new test is picked up by being
added — nothing to register.

## Main suite groups (`vmdhole/tests/run_tests.sh`)

Twenty groups. Each wrapper is a `test_<name>.sh`; the runner streams output,
carries each group's real exit status out of the pipeline, names a failing
group (`>>> <group>: FAILED (exit N)`), and lists skipped groups at the end.

| Group | Verifies |
|---|---|
| `test_headless_smoke` | the plugin loads, parses and imports under `vmd -dispdev text` — asserts numbers, not "it ran" |
| `test_accel_parity` | the shipped binaries ARE accelerated (OpenMP link + real scaling) and byte-identical to stock |
| `test_hydro_qco_parity` | the hydration C projection is bit-identical to the Tcl loop |
| `test_hole_tcl_fallback` | the pure-Tcl HOLE engine: same profile table, byte-identical `.sph`, refusals for cards it cannot honour |
| `test_hole_tcl_pore_methods` | the fallback under CONNOLLY and CAPSULE, including which table cells stay blank |
| `test_hole_tcl_fallback_e2e` | the fallback through `run_analysis` itself — per-frame file layout included |
| `test_hole_fast_coord` | the packed coordinate record is an accelerator, not a second answer |
| `test_capsule_incomplete` | HOLE's own early-stop warning reaches the user instead of being swallowed |
| `test_ellipse_parity` | the ellipse probe's C accelerator agrees with its Tcl reference |
| `test_h2dmap_parity` | the 2D-map parallel build: identical output and real speedup (the 1-thread cost is reported, not asserted) |
| `test_release_integrity` | the packaged ZIP is built fresh and self-consistent |
| `test_tcl_pitfalls` | repo lint (`tcl_lint.py`) plus Tcl-pitfall checks |
| `test_tunnel_separation` | pore and tunnel modes share no storage or output roots |
| `test_tunnel_clustering` | cross-frame cluster identity; display clustering leaves it untouched |
| `test_tunnel_import` | tunnel Save/Import round-trips byte-identically; combined HOLE+tunnel folders load |
| `test_mole_tcl_port` | the pure-Tcl MOLE engine reproduces the C engine slot for slot |
| `test_hcapen_cache` | HCAPEN's cutoff cache vs stock (auto-discovers `native/stock_build/hole2/src`, or pass a tree as `$1`) |
| `test_inline_current` | the inlined HOLE engine matches its `vmdhole/hole_tcl/` source |
| `test_adapter_schema` | the export adapter schema |
| `test_gui_reachable` | every Tunnel control reachable ON SCREEN in a real GUI, dialogs open/close/reopen, the lining window follows the selection, plus a second pass on a HET-carrying structure |

## Native verifications (`native/`)

| Check | What it proves |
|---|---|
| `verify.sh` Part A | the local `.sos` fixture corpus (`benchmarks/fixtures/`, gitignored; Part A skips when absent): stock vs fast triangulation byte-identical |
| `verify.sh` Part B | HOLE's own example structures through the full `hole -> sph_process -> sos_triangle` pipeline (the radius file is staged into the workdir — HOLE truncates long card paths silently) |
| `verify.sh` Part C | compiled `--hydro` colouring vs the Tcl reference (needs local, gitignored frames) |
| `verify.sh` Part D | `--points` dot surfaces: vertex sets identical |
| `verify.sh` Part E | parallel CONNOLLY / fast CAPSULE vs stock `hole`: `.sph`, profile and conductance byte-identical |
| `connolly_patches/test_hcapen_cache.sh` | the two known HCAPEN cache defects against a stock source tree |
| `hydration/test_hydro_project_unit.py` | `hydro_project.c` vs an independent Python reference, KDE support handling included |
| `hydration/test_hydro_accel_parity.tcl` | end-to-end hydration parity on real trajectories — its fixtures are not in the repository, so it must be run on a machine that has them |
