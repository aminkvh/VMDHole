#!/bin/sh
# Ion & Water > Openings: the scan samples every trajectory frame while HOLE
# may have analysed a subset. A sample on an unanalysed frame uses the nearest
# analysed frame's openings, so a stay across frames 0-5 with only frames 0
# and 4 analysed is ONE six-frame visit, not two one-frame visits.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "openings-occupancy"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

OUT=$(tclsh <<TCLEOF 2>&1
namespace eval ::VMDPathFinder {
    variable state; variable results; variable result_frames {0 4}
    variable plot_data_version 1; variable ion_flow_raw
    set results [dict create 0 [dict create sph_file "$TCL"] 4 [dict create sph_file "$TCL"]]
    proc _conn_site_table {} { return [dict create status ok frames {0 4} sites {{10 0 2 {}}}] }
    proc _conn_lobe_tol {} { return {6.0 35.0} }
    proc _conn_pore_margin {} { return 2.0 }
    proc _conn_frame_axis {fr} { return [list {0 0 1} {0 0 0} {}] }
    proc _conn_classify_cached {sph cv cp margin basis} { return [dict create lateral {{5.0 0.0 10.0}}] }
    proc _conn_frame_lobes {cls} { return {{0 0 0 {0}}} }
    proc _conn_lobe_site_map {table fr lobes} { return [dict create 0 7] }
    proc _conn_line_xyz {line} { return [lrange \$line 0 2] }
    proc _axg_unit {v} { return \$v }
    proc _conn_axis_basis {ux uy uz} { return {1 0 0 0 1 0} }
    proc _ion_flow_shell_value {} { return 3.0 }
}
$(awk '/^proc ::VMDPathFinder::_is_finite /,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_lobe_grid/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_cell_key/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_sample_site/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_nearest_int_in_sorted_list/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_open_close_episode/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_opening_occupancy/,/^}/' "$TCL")
namespace eval ::VMDPathFinder {
    # one Na+ in the opening on frames 0-5 (r 5 -> 8 A), one water on frames 2-3,
    # one Cl- far out in bulk (r 20 A) that must not count
    set ion_flow_raw [dict create has_water 1 nframes 6 traces [list \\
        [dict create species Na+ az {0 0 0 0 0 0} z {10 10 10 10 10 10} frame {0 1 2 3 4 5} r {5 5 5 5 5 8}] \\
        [dict create species Water az {0 0} z {10 10} frame {2 3} r {5 5}] \\
        [dict create species Cl- az {0 0} z {10 10} frame {0 1} r {20 20}]]]
    set occ [_conn_opening_occupancy]
    set e [dict get \$occ 7]
    puts "IONS [dict get \$e ions] VISITS [dict get \$e dwelln] DWELL [dict get \$e dwell] MOVED [dict get \$e cross] WATERS [dict get \$e waters] WDWELL [dict get \$e wdwell] SPECIES [dict get \$e species] SITES [dict size \$occ]"
}
TCLEOF
)
get() { printf '%s\n' "$OUT" | awk -v k="$1" '{ for (i = 1; i < NF; i++) if ($i == k) { print $(i+1); exit } }'; }
[ "$(get IONS)" = "6" ] && ok "six ion-frames counted across analysed and unanalysed frames" || bad "ions: $OUT"
[ "$(get VISITS)" = "1" ] && ok "...as one visit" || bad "visits: $(get VISITS)"
[ "$(get DWELL)" = "6.0" ] && ok "...six frames long" || bad "dwell: $(get DWELL)"
[ "$(get MOVED)" = "1" ] && ok "a 3 A radial change counts as moved" || bad "moved: $(get MOVED)"
[ "$(get WATERS)" = "2" ] && ok "water counted separately" || bad "waters: $(get WATERS)"
[ "$(get WDWELL)" = "2.0" ] && ok "...with its own visit length" || bad "wdwell: $(get WDWELL)"
printf '%s' "$OUT" | grep -q "SPECIES Na+ 6" && ok "species carry ion-frame counts" || bad "species: $OUT"
printf '%s' "$OUT" | grep -q "SITES 1" && ok "an ion in bulk beyond the opening is not attributed" || bad "sites: $OUT"
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
