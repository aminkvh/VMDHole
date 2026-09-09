#!/bin/sh
# Memory slots: several spherical analyses held at once, none overwriting another.
#
# The properties that matter, and that the UI will rest on:
#   * a new memory starts from DEFAULTS, not from a copy of the current one -
#     a second memory describes a DIFFERENT pore, so inheriting the last CPOINT
#     would be the wrong start;
#   * switching back restores that memory's own numbers AND its results;
#   * a new memory does not touch an older memory's results or work dir - the
#     whole point is that the first pore keeps describing the first pore;
#   * each memory runs into its own directory, so a second run cannot write
#     over the first one's frames;
#   * sync copies DISPLAY settings only, never run geometry.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "memory-slots"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

OUT=$(tclsh <<TCLEOF 2>&1
# The memory core is pure state handling - no VMD, no Tk. Pull in just those
# procs plus the variables they use.
namespace eval ::VMDPathFinder {
    variable state; variable default_state
    variable results [dict create]; variable result_frames {}
    variable pore_memories [dict create]; variable pore_memory_active ""
    variable pore_memory_next 1; variable plot_data_version 0
    variable pore_memory_runsig ""; variable last_geom_key ""
    variable mem_surface_mols; array set mem_surface_mols {}
    proc analysis_mode {} { return "pore" }
}
$(awk '/^proc ::VMDPathFinder::_mem_run_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_shared_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_stale_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_shared_signature/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_mark_fresh/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_stale /,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_stale_ids/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_track_ids/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_workdir_for/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_forget_tracks/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_display_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_enabled/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_capture/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_apply/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_stash_active/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_ensure_first/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_new/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_activate/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_delete/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_sync_presentation/,/^}/' "$TCL")

namespace eval ::VMDPathFinder {
    array set default_state {selection protein cpoint {} cvect {0 0 1} sample 0.25
                             endrad 15.0 surface_color hole_def pore_method circular
                             display_mode triangulated mesher csg search_engine mc
                             show_mean_surface 0 work_dir /tmp/base}
    array set state [array get default_state]
    set state(work_dir) /tmp/base
    # memory 1: a real pore
    set state(cpoint) "1 2 3"
    set state(selection) "protein and chain A"
    set state(surface_color) green
    set results [dict create 0 {a b}]
    set result_frames {0}
    _mem_ensure_first; _mem_mark_fresh
    set id2 [_mem_new]
    puts "ID2 \$id2"
    puts "NEW_CPOINT '\$state(cpoint)'"
    puts "NEW_SEL '\$state(selection)'"
    puts "NEW_RESULTS [dict size \$results]"
    puts "NEW_WORKDIR [_mem_workdir_for \$id2 \$state(work_dir)]"
    # fill memory 2 with its own pore
    set state(cpoint) "9 9 9"
    set results [dict create 5 {c d}]
    set result_frames {5}
    _mem_mark_fresh
    # back to 1
    _mem_activate 1
    puts "BACK_CPOINT '\$state(cpoint)'"
    puts "BACK_SEL '\$state(selection)'"
    puts "BACK_RESULTS [dict size \$results] frames \$result_frames"
    puts "BACK_COLOR \$state(surface_color)"
    # memory 2 must be untouched by having visited 1
    _mem_activate \$id2
    puts "M2_CPOINT '\$state(cpoint)'"
    puts "M2_RESULTS [dict size \$results] frames \$result_frames"
    # SHARED settings: a memory must not carry its own copy
    set _ov {}
    foreach _k [_mem_shared_keys] {
        if {[dict exists [dict get [dict get \$pore_memories 1] params] \$_k]} { lappend _ov \$_k }
    }
    puts "SHARED_IN_PARAMS [llength \$_ov] \$_ov"
    # sync presentation: display only
    set state(surface_color) red
    set state(dot_density) 30
    set n [_mem_sync_presentation 0]
    _mem_activate 1
    puts "SYNC_N \$n"
    puts "SYNC_COLOR \$state(surface_color)"
    puts "SYNC_DISPLAY \$state(dot_density)"
    puts "SYNC_CPOINT_UNCHANGED '\$state(cpoint)'"
    # work dirs must differ
    puts "WD1 [_mem_workdir_for 1 /tmp/base]"
    puts "WD2 [_mem_workdir_for \$id2 /tmp/base]"
    # a shared setting changing makes stored results stale
    puts "STALE_BEFORE [_mem_stale_ids]"
    set state(sample) 1.0
    puts "STALE_AFTER [_mem_stale_ids]"
    set state(sample) 0.25
    puts "STALE_RESTORED [_mem_stale_ids]"
    # cannot delete the last one
    _mem_delete 1
    puts "AFTER_DEL [lsort -integer [dict keys \$pore_memories]]"
    _mem_delete \$id2
    puts "AFTER_DEL2 [lsort -integer [dict keys \$pore_memories]]"
}
TCLEOF
)
get() { printf '%s\n' "$OUT" | sed -n "s/^$1 //p" | head -1; }

[ "$(get ID2)" = "2" ] && ok "the second memory gets its own id" || bad "id2=$(get ID2)"
[ "$(get NEW_CPOINT)" = "''" ] && ok "a new memory starts from the DEFAULT cpoint, not a copy" || bad "new cpoint $(get NEW_CPOINT)"
[ "$(get NEW_SEL)" = "'protein'" ] && ok "...and the default selection" || bad "new selection $(get NEW_SEL)"
[ "$(get NEW_RESULTS)" = "0" ] && ok "...with no results of its own yet" || bad "new results $(get NEW_RESULTS)"
case "$(get NEW_WORKDIR)" in */mem_2) ok "...and reports its own work dir" ;; *) bad "new workdir $(get NEW_WORKDIR)" ;; esac

[ "$(get BACK_CPOINT)" = "'1 2 3'" ] && ok "switching back restores that memory's cpoint" || bad "back cpoint $(get BACK_CPOINT)"
[ "$(get BACK_SEL)" = "'protein and chain A'" ] && ok "...and its selection" || bad "back selection $(get BACK_SEL)"
printf '%s' "$(get BACK_RESULTS)" | grep -q "^1 frames 0$" && ok "...and its own results, not the other memory's" || bad "back results: $(get BACK_RESULTS)"

[ "$(get M2_CPOINT)" = "'9 9 9'" ] && ok "the other memory is unchanged by the visit" || bad "m2 cpoint $(get M2_CPOINT)"
printf '%s' "$(get M2_RESULTS)" | grep -q "^1 frames 5$" && ok "...including its results" || bad "m2 results: $(get M2_RESULTS)"

[ "$(get SYNC_N)" = "1" ] && ok "sync reached the other memory" || bad "sync n=$(get SYNC_N)"
[ "$(get SYNC_COLOR)" = "red" ] && ok "sync copies the display colour" || bad "sync colour $(get SYNC_COLOR)"
[ "$(get SYNC_DISPLAY)" = "30" ] && ok "...and the dot density" || bad "sync display $(get SYNC_DISPLAY)"
[ "$(get SYNC_CPOINT_UNCHANGED)" = "'1 2 3'" ] && ok "sync does NOT touch run geometry" || bad "sync changed cpoint to $(get SYNC_CPOINT_UNCHANGED)"

W1=$(get WD1); W2=$(get WD2)
# The binding guarantee is asserted in test_gui_smoke.sh (MEMROOT), which
# exercises resolve_output_root - the function the run path actually calls.
[ -n "$W1" ] && [ "$W1" != "$W2" ] && ok "the work-dir helper reports a different directory per memory ($W1 vs $W2)" || bad "work dirs collide: '$W1' '$W2'"

SIP=$(get SHARED_IN_PARAMS)
[ "${SIP%% *}" = "0" ] && ok "no shared setting is stored per memory (sampling/mesher/engine cannot be mixed)" || bad "shared keys leaked into a memory: $SIP"
[ "$(get STALE_BEFORE)" = "" ] && ok "results match the settings they were produced under" || bad "stale before: $(get STALE_BEFORE)"
printf '%s' "$(get STALE_AFTER)" | grep -q "1" && ok "changing the sampling marks stored results stale" || bad "stale after: '$(get STALE_AFTER)'"
[ "$(get STALE_RESTORED)" = "" ] && ok "...and putting it back clears the flag" || bad "stale restored: $(get STALE_RESTORED)"

[ "$(get AFTER_DEL)" = "2" ] && ok "a memory can be deleted" || bad "after delete: $(get AFTER_DEL)"
[ "$(get AFTER_DEL2)" = "2" ] && ok "the LAST memory cannot be deleted" || bad "after deleting the last: $(get AFTER_DEL2)"

echo "memory-slots: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || { echo "--- raw ---"; printf '%s\n' "$OUT" | head -25; exit 1; }
