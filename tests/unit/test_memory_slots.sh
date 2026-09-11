#!/bin/sh
# Memory slots: each memory is one run, held with its settings, its results and
# the folder it lives in.
#   * a new memory starts from the defaults for the pore's own settings, and
#     keeps the method settings and the frame range;
#   * switching back restores that memory's settings, results and folder;
#   * a new memory has no folder until it runs (the run makes one);
#   * sync copies display settings only, never run geometry;
#   * a memory whose run does not cover the shown frame is blanked, not left
#     showing the last frame it did cover;
#   * the saved settings (not the shipped ones) are what a new memory starts from.
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
namespace eval ::VMDPathFinder {
    proc _sync_point_marker {args} {}
    proc _mem_point_surface_mol {} {}
    proc _mem_seed_axis {} {}
    proc _mem_keep_frame {} {}
    proc _sync_cvect_handles {args} {}
    proc _mem_sync_shared_labels {} {}
    variable state; variable default_state
    variable results [dict create]; variable result_frames {}
    variable pore_memories [dict create]; variable pore_memory_active ""
    variable pore_memory_next 1; variable plot_data_version 0
    variable run_root ""; variable run_root_temp 0; variable last_geom_key ""
    variable current_surface_mol -1
    variable mem_surface_mols; array set mem_surface_mols {}
    proc analysis_mode {} { return "hole" }
}
$(awk '/^proc ::VMDPathFinder::_mem_run_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_track_ids/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_forget_tracks/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_sync_cues/,/^}/' "$TCL")
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
$(awk '/^proc ::VMDPathFinder::_config_persistent_keys/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_config_adopt_defaults/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_render_other_memories/,/^}/' "$TCL")

namespace eval ::VMDPathFinder {
    array set default_state {selection protein cpoint {} cvect {0 0 1} sample 0.25
                             endrad 15.0 surface_color hole_def pore_method circular
                             display_mode triangulated mesher csg search_engine mc
                             conn_engine hole show_mean_surface 0 frame_spec now}
    array set state [array get default_state]
    # the user's saved default endrad is 20, the shipped one 15
    set state(endrad) 20.0
    _config_adopt_defaults
    # memory 1: a real run
    set state(cpoint) "1 2 3"
    set state(selection) "protein and chain A"
    set state(surface_color) green
    set state(pore_method) connolly
    set state(frame_spec) all
    set results [dict create 0 {a b}]
    set result_frames {0}
    set run_root /tmp/base/1BL8_20260910-120000-abc123
    _mem_ensure_first
    set id2 [_mem_new]
    puts "ID2 \$id2"
    puts "NEW_CPOINT '\$state(cpoint)'"
    puts "NEW_SEL '\$state(selection)'"
    puts "NEW_METHOD \$state(pore_method) \$state(frame_spec)"
    puts "NEW_RESULTS [dict size \$results]"
    puts "NEW_ROOT '\$run_root'"
    puts "NEW_ENDRAD \$state(endrad)"
    # memory 2 runs its own pore into its own folder
    set state(cpoint) "9 9 9"
    set results [dict create 5 {c d}]
    set result_frames {5}
    set run_root /tmp/base/1BL8_20260910-130000-def456
    _mem_stash_active
    _mem_activate 1
    puts "BACK_CPOINT '\$state(cpoint)'"
    puts "BACK_SEL '\$state(selection)'"
    puts "BACK_RESULTS [dict size \$results] frames \$result_frames"
    puts "BACK_COLOR \$state(surface_color)"
    puts "BACK_ROOT '\$run_root'"
    _mem_activate \$id2
    puts "M2_CPOINT '\$state(cpoint)'"
    puts "M2_RESULTS [dict size \$results] frames \$result_frames"
    puts "M2_ROOT '\$run_root'"
    # sync presentation: display only
    set state(surface_color) red
    set state(dot_density) 30
    set n [_mem_sync_presentation 0]
    _mem_activate 1
    puts "SYNC_N \$n"
    puts "SYNC_COLOR \$state(surface_color)"
    puts "SYNC_DISPLAY \$state(dot_density)"
    puts "SYNC_CPOINT_UNCHANGED '\$state(cpoint)'"
    # cannot delete the last one
    _mem_delete 1
    puts "AFTER_DEL [lsort -integer [dict keys \$pore_memories]]"
    _mem_delete \$id2
    puts "AFTER_DEL2 [lsort -integer [dict keys \$pore_memories]]"
}

# Drawing the frame for every OTHER memory: one covers frame 3, one does not.
namespace eval ::VMDPathFinder {
    proc _abort_requested {} { return 0 }
    proc _op_in_progress {} { return 0 }
    proc load_surface_for_frame {frame {draft 0}} { lappend ::DREW "mem\$::VMDPathFinder::pore_memory_active f\$frame" }
    proc _mem_blank_track {frame} { lappend ::BLANKED "mem\$::VMDPathFinder::pore_memory_active f\$frame" }
    set pore_memories [dict create \
        1 [dict create params {surface_color red} results [dict create 0 {} 1 {} 2 {} 3 {}] frames {0 1 2 3} run_root /tmp/a] \
        2 [dict create params {surface_color blue} results [dict create 0 {} 1 {} 2 {}] frames {0 1 2} run_root /tmp/b] \
        3 [dict create params {surface_color green} results [dict create 0 {} 1 {} 2 {}] frames {0 1 2} run_root /tmp/c]]
    set pore_memory_active 3
    set state(pore_method) circular
    set ::DREW {}; set ::BLANKED {}
    _mem_render_other_memories 3 1
    puts "DREW3 \$::DREW"
    puts "BLANKED3 \$::BLANKED"
    puts "ACTIVE_AFTER \$pore_memory_active"
}
TCLEOF
)
get() { printf '%s\n' "$OUT" | sed -n "s/^$1 //p" | head -1; }

[ "$(get ID2)" = "2" ] && ok "the second memory gets its own id" || bad "id2=$(get ID2)"
[ "$(get NEW_CPOINT)" != "'1 2 3'" ] && ok "a new memory never inherits the previous memory's cpoint" || bad "new cpoint $(get NEW_CPOINT)"
[ "$(get NEW_SEL)" = "'protein'" ] && ok "...and starts from the default selection" || bad "new selection $(get NEW_SEL)"
[ "$(get NEW_METHOD)" = "connolly all" ] && ok "...but keeps the method and the frame range" || bad "new method/frames $(get NEW_METHOD)"
[ "$(get NEW_RESULTS)" = "0" ] && ok "...with no results of its own yet" || bad "new results $(get NEW_RESULTS)"
[ "$(get NEW_ROOT)" = "''" ] && ok "...and no folder until it runs" || bad "new root $(get NEW_ROOT)"
[ "$(get NEW_ENDRAD)" = "20.0" ] && ok "...from the saved defaults, not the shipped ones" || bad "new endrad $(get NEW_ENDRAD)"

[ "$(get BACK_CPOINT)" = "'1 2 3'" ] && ok "switching back restores that memory's cpoint" || bad "back cpoint $(get BACK_CPOINT)"
[ "$(get BACK_SEL)" = "'protein and chain A'" ] && ok "...and its selection" || bad "back selection $(get BACK_SEL)"
printf '%s' "$(get BACK_RESULTS)" | grep -q "^1 frames 0$" && ok "...and its own results" || bad "back results: $(get BACK_RESULTS)"
[ "$(get BACK_ROOT)" = "'/tmp/base/1BL8_20260910-120000-abc123'" ] && ok "...and its own folder" || bad "back root $(get BACK_ROOT)"

[ "$(get M2_CPOINT)" = "'9 9 9'" ] && ok "the other memory is unchanged by the visit" || bad "m2 cpoint $(get M2_CPOINT)"
printf '%s' "$(get M2_RESULTS)" | grep -q "^1 frames 5$" && ok "...including its results" || bad "m2 results: $(get M2_RESULTS)"
[ "$(get M2_ROOT)" = "'/tmp/base/1BL8_20260910-130000-def456'" ] && ok "...and its folder" || bad "m2 root $(get M2_ROOT)"

[ "$(get SYNC_N)" = "1" ] && ok "sync reached the other memory" || bad "sync n=$(get SYNC_N)"
[ "$(get SYNC_COLOR)" = "red" ] && ok "sync copies the display colour" || bad "sync colour $(get SYNC_COLOR)"
[ "$(get SYNC_DISPLAY)" = "30" ] && ok "...and the dot density" || bad "sync display $(get SYNC_DISPLAY)"
[ "$(get SYNC_CPOINT_UNCHANGED)" = "'1 2 3'" ] && ok "sync does not touch run geometry" || bad "sync changed cpoint to $(get SYNC_CPOINT_UNCHANGED)"

[ "$(get AFTER_DEL)" = "2" ] && ok "a memory can be deleted" || bad "after delete: $(get AFTER_DEL)"
[ "$(get AFTER_DEL2)" = "2" ] && ok "the last memory cannot be deleted" || bad "after deleting the last: $(get AFTER_DEL2)"

[ "$(get DREW3)" = "{mem1 f3}" ] && ok "a memory whose run covers the frame is drawn" || bad "drew: $(get DREW3)"
[ "$(get BLANKED3)" = "{mem2 f3}" ] && ok "...and one whose run does not is blanked, not left stale" || bad "blanked: $(get BLANKED3)"
[ "$(get ACTIVE_AFTER)" = "3" ] && ok "...with the active memory restored either way" || bad "active after: $(get ACTIVE_AFTER)"

echo "memory-slots: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || { echo "--- raw ---"; printf '%s\n' "$OUT" | head -25; exit 1; }
