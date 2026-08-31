#!/bin/sh
# Save/Import for tunnel mode (Phase 2): reconstructing tunnel_results from
# each frame's out.dat, and finding it on disk in the first place.
#
#   1 round-trip: clobber live tunnel_results/tunnel_lining, reimport from
#     out.dat, and the reconstruction must be byte-identical to what the
#     live run actually computed - not merely non-empty.
#   2 a frame whose out.dat is missing is SKIPPED and named in the status
#     message, never silently absent without a trace and never backfilled
#     from .sph (whose centreline is trimmed differently for display -
#).
#   3 import_all_results_from_folder finds tunnel data with zero HOLE data
#     present and does not error or invent a spurious HOLE-side failure.
#   4 cluster-keyed state (tunnel_shown_cid, tunnel_gear_cid, the pinned
#     selection, the shown-default) is reset by import, not just its
#     rank-keyed mirrors - cids are renumbered from scratch on every import
#     (sequential, same as rank), so a stale entry would otherwise silently
#     reattach to whatever route now gets that number, and a stale "hide by
#     default" would hide every freshly imported route outright.
#   5 the same cluster-keyed gear mirror is reset by a RE-RUN of
#     run_tunnel_analysis too (check 4's twin - the same defect existed in
#     both the run and import reset paths).
#   6 the unified import entry point finds BOTH HOLE and tunnel results when
#     a real run of each shares one common work_dir (resolve_output_root's
#     own nested-not-sibling layout for an explicit work_dir).

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SRC="$DIR/../vmdhole.tcl"
ENGINE="$DIR/../../native/mole_tunnel_engine"
PDB="$DIR/fixtures/mole_reference/1BL8.pdb"
echo "=============================================================="
echo "tunnel-import: $SRC"
pass=0; fail=0
ok()  { echo "  PASS  $1"; pass=$((pass+1)); }
bad() { echo "  FAIL  $1"; fail=$((fail+1)); }

[ -f "$SRC" ] || { echo "  FAIL  source not found"; exit 1; }
if [ ! -x "$ENGINE" ] || [ ! -f "$PDB" ]; then
    # Group-level ^SKIP: form - see test_tunnel_clustering.sh for why.
    echo "SKIP: no compiled mole_tunnel_engine or no 1BL8 fixture"
    exit 0
fi

VMD=${VMD:-vmd}
if ! command -v "$VMD" >/dev/null 2>&1; then
    # Nothing below runs without vmd; the silenced invocation used to leave no
    # output file and the group failed with "script produced no output" on any
    # machine without vmd, instead of skipping like every other vmd group.
    echo "SKIP: no vmd on PATH"
    exit 0
fi
OUT=$(mktemp)
WORK=$(mktemp -d)
cat > "$OUT.tcl" <<EOF
package provide Tk 8.5
if {[catch {source "$SRC"} e]} {
    set o [open "$OUT.out" w]; puts \$o "FAIL source: \$e"; puts \$o "BAD=1"; close \$o; quit
}
set o [open "$OUT.out" w]
set bad 0

set mid [mol new "$PDB" waitfor all]
animate dup \$mid
animate dup \$mid

set ::VMDHole::state(molid) \$mid
set ::VMDHole::state(selection) "protein"
set ::VMDHole::state(engine) "mole"
set ::VMDHole::state(tunnel_start) "73.853 26.536 26.594"
set ::VMDHole::state(mole_engine_exec) "$ENGINE"
set ::VMDHole::state(work_dir) "$WORK"
set ::VMDHole::state(frame_spec) "all"
::VMDHole::run_tunnel_analysis
set nfr [llength \$::VMDHole::tunnel_result_frames]
if {\$nfr < 2} {
    puts \$o "FAIL run produced \$nfr frame(s), need >=2: \$::VMDHole::state(status)"
    incr bad
} else {
    set root \$::VMDHole::tunnel_root
    set orig_frames \$::VMDHole::tunnel_result_frames
    array set orig_results {}
    foreach fr \$orig_frames { set orig_results(\$fr) \$::VMDHole::tunnel_results(\$fr) }

    # ---- 1. round-trip byte-identical ----
    array unset ::VMDHole::tunnel_results; array set ::VMDHole::tunnel_results {}
    set ::VMDHole::tunnel_result_frames {}
    ::VMDHole::import_tunnel_results_from_folder \$root
    set match 1
    if {\$::VMDHole::tunnel_result_frames ne \$orig_frames} { set match 0 }
    foreach fr \$orig_frames {
        if {![info exists ::VMDHole::tunnel_results(\$fr)] \
                || \$::VMDHole::tunnel_results(\$fr) ne \$orig_results(\$fr)} { set match 0 }
    }
    if {\$match} { puts \$o "PASS round-trip byte-identical" } else {
        puts \$o "FAIL round-trip mismatch (frames=\$::VMDHole::tunnel_result_frames want=\$orig_frames)"
        incr bad
    }

    # ---- 2. missing out.dat is skipped and named ----
    set fr1 [lindex \$orig_frames 1]
    set fd1 [file join \$root [format "tunnel_%05d" \$fr1]]
    file rename [file join \$fd1 out.dat] [file join \$fd1 out.dat.bak]
    array unset ::VMDHole::tunnel_results; array set ::VMDHole::tunnel_results {}
    set ::VMDHole::tunnel_result_frames {}
    ::VMDHole::import_tunnel_results_from_folder \$root
    set skip_ok [expr {![info exists ::VMDHole::tunnel_results(\$fr1)] \
        && [lsearch -exact \$::VMDHole::tunnel_result_frames \$fr1] < 0 \
        && [string match "*skipped*\$fr1*" \$::VMDHole::state(status)]}]
    if {\$skip_ok} { puts \$o "PASS missing out.dat skipped and named in status" } else {
        puts \$o "FAIL missing out.dat not handled correctly: \$::VMDHole::state(status)"
        incr bad
    }
    file rename [file join \$fd1 out.dat.bak] [file join \$fd1 out.dat]

    # ---- 3. unified import on a tunnel-only directory ----
    set ::VMDHole::state(import_dir) \$root
    array unset ::VMDHole::tunnel_results; array set ::VMDHole::tunnel_results {}
    set ::VMDHole::tunnel_result_frames {}
    if {[catch {::VMDHole::import_all_results_from_folder} uerr]} {
        puts \$o "FAIL import_all_results_from_folder errored on tunnel-only dir: \$uerr"
        incr bad
    } elseif {[llength \$::VMDHole::tunnel_result_frames] == 0} {
        puts \$o "FAIL import_all_results_from_folder found no tunnel data"
        incr bad
    } else {
        puts \$o "PASS unified import finds tunnel-only data cleanly"
    }

    # ---- 4. cluster-keyed state is reset by import, not just the rank-keyed
    # arrays (see the header comment above) ----
    set ::VMDHole::tunnel_shown_cid(1) 0
    set ::VMDHole::tunnel_gear_cid(1,color) red
    set ::VMDHole::state(tunnel_selected_cid) 999
    set ::VMDHole::state(tunnel_shown_default) 0
    array unset ::VMDHole::tunnel_results; array set ::VMDHole::tunnel_results {}
    set ::VMDHole::tunnel_result_frames {}
    ::VMDHole::import_tunnel_results_from_folder \$root
    set cleared [expr {![info exists ::VMDHole::tunnel_shown_cid(1)] \
        && ![info exists ::VMDHole::tunnel_gear_cid(1,color)] \
        && \$::VMDHole::state(tunnel_selected_cid) eq "" \
        && \$::VMDHole::state(tunnel_shown_default) == 1}]
    if {\$cleared} { puts \$o "PASS import clears cluster-keyed shown/gear state and the stale selection" } else {
        puts \$o "FAIL import left stale cluster-keyed state (shown_cid/gear_cid/selected_cid/shown_default not all reset)"
        incr bad
    }

    # ---- 5. the same cluster-keyed gear mirror is reset by a RE-RUN of
    # run_tunnel_analysis too, not just by import (check 4's twin) ----
    set ::VMDHole::tunnel_gear_cid(1,color) red
    ::VMDHole::run_tunnel_analysis
    if {[info exists ::VMDHole::tunnel_gear_cid(1,color)]} {
        puts \$o "FAIL re-run left a stale tunnel_gear_cid(1,color) from the previous run"
        incr bad
    } else {
        puts \$o "PASS re-run of run_tunnel_analysis clears the cluster-keyed gear mirror"
    }
}

# ---- 4. combined HOLE + tunnel from one common work_dir ----
set WORK2 "${WORK}_combo"
file mkdir \$WORK2
set ::VMDHole::state(cpoint) "73.853 26.536 26.594"
set ::VMDHole::state(work_dir) \$WORK2
set ::VMDHole::state(frame_spec) "now"
::VMDHole::run_analysis
::VMDHole::run_tunnel_analysis
set ::VMDHole::results [dict create]
set ::VMDHole::result_frames {}
array unset ::VMDHole::tunnel_results; array set ::VMDHole::tunnel_results {}
set ::VMDHole::tunnel_result_frames {}
set ::VMDHole::state(import_dir) \$WORK2
::VMDHole::import_all_results_from_folder
set both [expr {[llength \$::VMDHole::result_frames] > 0 && [llength \$::VMDHole::tunnel_result_frames] > 0}]
if {\$both} { puts \$o "PASS combined HOLE+tunnel import from one common directory" } else {
    puts \$o "FAIL combined import: HOLE frames=\$::VMDHole::result_frames tunnel frames=\$::VMDHole::tunnel_result_frames"
    incr bad
}

puts \$o "BAD=\$bad"
close \$o
quit
EOF
"$VMD" -dispdev text -e "$OUT.tcl" >/dev/null 2>&1
if [ -f "$OUT.out" ]; then
    while IFS= read -r line; do
        case "$line" in
            PASS*) ok "${line#PASS }" ;;
            FAIL*) bad "${line#FAIL }" ;;
        esac
    done < "$OUT.out"
    if ! grep -q "^BAD=" "$OUT.out"; then
        bad "script did not complete (no BAD= sentinel)"
    fi
else
    bad "script produced no output"
fi
rm -rf "$OUT" "$OUT.tcl" "$OUT.out" "$WORK" "${WORK}_combo"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
