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
| `test_skip_message_contract.sh` | Group-level skips written in a form `run_tests.sh`'s `^SKIP:` gate cannot see, which makes a group that checked nothing count as a PASS — and pass the `VMDHOLE_RELEASE=1` gate. |

## Wiring these into the main suite

They are intentionally *not* wired in yet. To adopt one, add its name to
`GROUPS` in `vmdhole/tests/run_tests.sh` and give it a `test_*.sh` wrapper in
`vmdhole/tests/`, or call `run_unit_tests.sh` as a single group. Both of these
are also good candidates for the `tcl-lint` / `build-c` CI jobs in
`.github/workflows/tests.yml`, because neither needs VMD.
