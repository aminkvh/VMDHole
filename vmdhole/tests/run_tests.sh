#!/bin/sh
# VMDHole regression suite.
#
# Each test targets a defect class that ACTUALLY occurred in this project and was
# found by hand rather than by a check. Every one has been verified to go red on
# the real defect before being trusted - a green suite that could not have caught
# the bugs it is named after is worse than no suite.
#
#   test_headless_smoke     the plugin must load, parse and import under
#                           `vmd -dispdev text`, which has NO Tk. Asserts numbers,
#                           not "it ran". Catches: unguarded winfo/tk_messageBox.
#   test_accel_parity       the shipped binaries must BE accelerated (Makefile
#                           invariants + OpenMP link + real scaling) AND still be
#                           byte-identical to stock. Catches: the R-15 class, where
#                           a build silently produced a SERIAL sph_process with
#                           perfectly correct output.
#   test_release_integrity  builds packaging/vmdhole.zip fresh from git ls-files
#                           and checks it: no file shipped twice, pkgIndex
#                           points somewhere real, version self-consistent.
#   test_hole_tcl_fallback  the pure-Tcl HOLE engine, used when no `hole` binary
#                           is installed, must produce the SAME profile table and
#                           a byte-identical .sph - and refuse any control card it
#                           cannot honour. Catches: a fallback that silently drops
#                           a card (ignore, an unknown card) and returns a
#                           plausible profile for a run the user did not ask for.
#   test_hole_tcl_pore_methods
#                           the same fallback under CONNOLLY and CAPSULE, whose
#                           profile tables are not the spherical one with other
#                           numbers (CAPSULE is hcapgr.f; CONNOLLY adds three
#                           columns that must stay BLANK on mid-point rows).
#                           Needs the reference binary; skips without it.
#   test_hole_fast_coord    the packed coordinate record must be an accelerator,
#                           not a second answer: same .sph as the PDB round trip,
#                           and stock HOLE / the Tcl engine still get a PDB.
#   test_hole_tcl_fallback_e2e
#                           the same fallback through run_analysis itself: the
#                           pieces working is not the same as the plugin running,
#                           because the fallback has to reproduce the per-frame
#                           FILE layout the job pool's parser reads.
#   test_capsule_incomplete HOLE can stop its axial search early, SAY SO, and still
#                           exit 0 - so the exit-status check cannot see it and the
#                           frame becomes a normal result with a short profile.
#                           Catches: the plugin swallowing the engine's own warning,
#                           which is what "the capsule lining is far from the pore"
#                           was. Carries its own spherical control, so a detector
#                           that fired on every run would fail it.
#   test_ellipse_parity     the ellipse probe's C accelerator must agree with its
#                           Tcl reference. The LAST of the four engine pairs to
#                           get coverage - test_accel_parity compares C against
#                           STOCK C, never against Tcl. Catches: an accelerator
#                           that diverges from the reference it replaced.
#   test_mole_tcl_port      the pure-Tcl MOLE engine must reproduce the C engine
#                           slot for slot, not merely produce a valid Delaunay
#                           triangulation. Catches: a stage that looks like a
#                           transcription and is not, whose error surfaces
#                           several stages later as a plausible wrong tunnel.
#   test_gui_reachable      every Tunnel control must be reachable ON SCREEN,
#                           not merely created. SKIPS without a display.
#                           Catches: a panel taller than its sidebar, where
#                           `winfo ismapped` is 1 for controls the user cannot
#                           scroll to - which is what happened to the custom
#                           exit and path fields.
#
#   test_h2dmap_parity      the PARALLEL 2DMAPS routine must be byte-identical to
#                           the serial one at any thread count, AND actually be
#                           parallel. Catches both halves: a directive that
#                           changes the arithmetic, and a missing per-file
#                           -fopenmp rule, which leaves the directives as
#                           comments and builds a correct but serial routine.
#                           Needs a hole2/src checkout; skips without one.
#   test_inline_current     the HOLE engine inlined in vmdhole.tcl is GENERATED
#                           from vmdhole/hole_tcl/ - re-inline into a copy and require
#                           byte-identity. Catches: editing one copy and not the
#                           other, after which "the source of truth" is whichever
#                           file the reader happens to open.
#   test_hcapen_cache       the accelerated CAPSULE routine's cutoff cache must
#                           agree with stock HOLE across CALLS - both its known
#                           defects (a truncated candidate list, and reuse of a
#                           previous dataset's geometry) live in SAVE'd state,
#                           so no single-shot fixture can see them. verify.sh
#                           missed both. Needs a hole2/src checkout.
#
# Usage:  ./run_tests.sh          (exit 0 = all passed)
#         VMDHOLE_RELEASE=1 ./run_tests.sh
#                           release gate: a test that SKIPPED because its input
#                           was absent counts as a FAILURE. Several scientific
#                           oracles (real trajectories, a hole2/src checkout,
#                           CHAP output) are gitignored, so on a fresh clone
#                           those groups skip and the suite still exits 0 -
#                           green, having checked nothing. Do not tag a release
#                           from a run that is not green in this mode.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
fails=0; ran=0; skipped=""; missing=""
RELEASE="${VMDHOLE_RELEASE:-0}"
GROUPS="test_headless_smoke test_accel_parity test_hydro_qco_parity test_hole_tcl_fallback test_hole_tcl_pore_methods test_hole_tcl_fallback_e2e test_hole_fast_coord test_capsule_incomplete test_ellipse_parity test_h2dmap_parity test_release_integrity test_tcl_pitfalls test_tunnel_separation test_tunnel_clustering test_tunnel_import test_mole_tcl_port test_hcapen_cache test_inline_current test_adapter_schema test_gui_reachable"
EXPECTED=$(echo $GROUPS | wc -w)
for t in $GROUPS; do
    # A group whose file is gone or lost its exec bit used to vanish silently -
    # the suite still printed ALL <fewer> GROUPS PASSED. Name it and fail.
    if [ ! -x "$DIR/$t.sh" ]; then missing="$missing $t"; continue; fi
    ran=$((ran+1))
    echo "=============================================================="
    # Streams as it runs AND captures, so a 20-minute group still shows
    # progress. The __RC__ marker carries the real exit status back out of the
    # pipeline (POSIX sh has no PIPESTATUS) and is filtered from the display.
    _log=$(mktemp)
    ( "$DIR/$t.sh" 2>&1; echo "__RC__$?" ) | tee "$_log" | grep -v '^__RC__'
    _rc=$(sed -n 's/^__RC__//p' "$_log")
    [ "${_rc:-1}" -eq 0 ] || fails=$((fails+1))
    # Anchored to the GROUP-level form ("SKIP: no vmd on PATH"), not a bare
    # substring. An unanchored match also caught the per-assertion form
    # ("  SKIP  Ion Flow's axial window - this structure has no ions"), which
    # promoted a group that ran in full to "checked nothing" - and under
    # VMDHOLE_RELEASE=1 that is a FAILURE, so one benign internal skip could
    # block a release. test_gui_reachable and test_headless_smoke were both
    # mislabelled this way.
    if grep -q "^SKIP:" "$_log"; then skipped="$skipped $t"; fi
    rm -f "$_log"
done
echo "=============================================================="
if [ -n "$missing" ]; then
    echo "MISSING or not executable:$missing"
    echo "  ^ expected $EXPECTED groups, ran $ran. A group that does not run is not a pass."
    fails=$((fails+1))
fi
if [ -n "$skipped" ]; then
    echo "SKIPPED (input absent):$skipped"
    if [ "$RELEASE" != "0" ]; then
        echo "VMDHOLE_RELEASE set - a skipped group is a FAILURE for a release build."
        fails=$((fails+1))
    else
        echo "  ^ these checked nothing. Re-run with VMDHOLE_RELEASE=1 before tagging."
    fi
fi
if [ "$fails" -eq 0 ] && [ "$ran" -eq "$EXPECTED" ]; then echo "ALL $ran TEST GROUPS PASSED"; else echo "$fails of $ran TEST GROUPS FAILED"; fi
[ "$fails" -eq 0 ]
