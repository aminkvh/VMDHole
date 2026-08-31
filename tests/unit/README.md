# `tests/unit/` — small, self-contained regression tests

Deliberately separate from the two existing suites so they can be wired into the
right runner on purpose rather than by accident:

| suite | scope | needs |
|---|---|---|
| `vmdhole/tests/` | the plugin, end to end | VMD, locally built HOLE binaries, fixtures |
| `vmdhole/hole_tcl/tests/` | the pure-Tcl HOLE engine vs a reference binary | reference `hole` build |
| **`tests/unit/`** | **one defect each, in isolation** | **a C compiler and/or `sh` — nothing else** |

## The rules these follow

1. **No heavy prerequisites.** No VMD, no trajectory data, no gitignored fixture
   corpus, no network. Anything needed is generated into a temp dir and cleaned
   up. That is what makes them safe to run on every commit, unlike the main
   suite whose heavier groups legitimately skip on a bare machine.
2. **Each is verified to go RED on the real defect before being trusted.** A
   green test that could not have caught the bug it is named after is worse
   than no test. Every file below records the exact defect it guards and the
   sanitizer output or symptom it was proven against.
3. **They test the contract, not the incident.** Where a bug is an instance of a
   class, the test enforces the class. `test_skip_message_contract.sh` was
   written for one file and immediately found three more instances.
4. **A test that genuinely cannot run prints `SKIP:` at column 0 and exits 0** —
   the same contract `vmdhole/tests/run_tests.sh` uses.

## Running

```sh
sh tests/unit/run_unit_tests.sh      # all of them; exit 0 = all passed
sh tests/unit/test_sos_nonfinite_input.sh    # or one at a time
```

## What is here

| test | guards |
|---|---|
| `test_sos_nonfinite_input.sh` | Signed-integer overflow (UB) in `sos_triangle_fast`'s two spatial grids when a `.sos` carries a non-finite coordinate. Introduced by the fork's spatial-hash optimisation — upstream has no cell index at all. Also asserts the fix did not change valid-input output, since byte-identity is this codebase's whole justification. |
| `test_sos_short_record_normals.sh` | A short `.sos` point record inheriting the PREVIOUS record's normal (or an indeterminate stack value on the first record), which `reorder_triangle` then feeds into vertex-order decisions. A short record must be indistinguishable from an explicit `0 0 0` normal. |
| `test_skip_message_contract.sh` | Group-level skips written in a form `run_tests.sh`'s `^SKIP:` gate cannot see, which makes a group that checked nothing count as a PASS — and pass the `VMDHOLE_RELEASE=1` gate. |
| `test_mole_cli_validation.sh` | The tunnel engine's CLI crashing on or silently reinterpreting nonsense: a negative `--cover` recursing to a heap overflow, an out-of-range weight running the default, a misplaced flag parsed as the probe radius. |
| `test_mole_coord_validation.sh` | `nan`/`inf`/junk atom coordinates parsing as 0.0 (or propagating into the predicates) and producing a plausible tunnel set with exit 0 instead of a refusal. |
| `test_mole_depth_bound.sh` | `remove_shallow` iterating the raw unclamped MinDepth value — linear cost in the parameter, ~12 h at `INT_MAX` — instead of the depth levels that occur in the graph. |
| `test_preset_exec_paths.sh` | `load_config` overwriting engine paths a batch script set before `init_executables`, silently selecting the ~100x slower Tcl fallback. |
| `test_headless_scratch_reclaim.sh` | Scratch-directory reclaim being reachable only through the GUI, so headless jobs killed by a scheduler leaked `/dev/shm` directories until reboot. |
| `test_headless_run_guards.sh` | Three headless run-path guards: a non-numeric or shell-metacharacter start point must be refused before it reaches `run.sh`/`exec sh -c`; an error escaping `run_tunnel_analysis` must restore `busy` and balance `_end_calc`; the empty-selection branch in `run_analysis` must not delete the atomselect handle shared across frames. |
| `test_tsv_reader_parity.sh` | The threaded TSV reader skipping `_resolve_conn_radii`, so a CONNOLLY trajectory's bottleneck differed 2.5x between the serial and threaded paths. |
| `test_gui_smoke.sh` | The GUI itself, without VMD: vmdhole.tcl sourced under plain tclsh+Tk with VMD stubbed, the real widget tree built, and 24 scripted user actions asserting the close-path (abort survives an in-flight close, traces survive reopen, a pending modal is answered), the busy guards, the deleted-molecule dialog class, the tunnel gear-popup route pinning, and nan-tolerant option fields. Skips without Tk or a display (CI runs it under xvfb). |
| `test_tsv_publish_on_failure.sh` | A failed profile parse (zero rows, or a throw) truncating the good `hole_profile.tsv` beside it. Both TSV writers must publish by rename only after a parse that produced rows. |

## Wiring

CI runs the whole directory: the `tcl-suite` job in
`.github/workflows/tests.yml` calls `run_unit_tests.sh` (it needs only tclsh,
a C compiler and python3, all of which that job installs). The main suite does
not invoke it; to adopt a test as a `run_tests.sh` group, add its name to
`GROUPS` there and give it a `test_*.sh` wrapper in `vmdhole/tests/`.
