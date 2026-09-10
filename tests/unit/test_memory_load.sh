#!/bin/sh
# Loading memories: every run folder under a chosen folder becomes a memory,
# in name order, and a run made on another structure is told apart by the
# atom count and radius of gyration its signature recorded.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "memory-load"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
mksig() { # dir natoms rgyr
  mkdir -p "$1/frame_0"
  printf 'sel|protein|cpoint|1 2 3|sample|0.25|endrad|15.0|method|circular|source|1BL8.pdb|coords|%s:1:72.187,26.600,19.890:%s|search|mc|conneng|hole|plugin|1.0' \
    "$2" "$3" > "$1/frame_0/vmdpathfinder_signature.dat"
}
# one structure's output folder holding three runs
mksig "$TMP/hole_output_1BL8/1BL8_20260910-120000-aaaaaa" 2824 20.607
mksig "$TMP/hole_output_1BL8/1BL8_20260910-130000-bbbbbb" 2824 20.607
mksig "$TMP/hole_output_1BL8/1BL8_20260910-140000-cccccc" 2824 20.607
mkdir -p "$TMP/hole_output_1BL8/notes"
# a plain run folder chosen directly
mksig "$TMP/plain" 2824 20.607
# a run made on a different structure
mksig "$TMP/other" 4649 25.100

OUT=$(tclsh <<TCLEOF 2>&1
proc lmapish {pairs} { set o {}; foreach p \$pairs { lappend o [lindex \$p 0] }; return \$o }
proc dirs {pairs} { set o {}; foreach p \$pairs { lappend o [file tail [lindex \$p 1]] }; return \$o }
# a stand-in molecule: 2824 atoms, rgyr 20.607 (frame 0)
proc molinfo {id what args} { if {\$what eq "get"} { switch -- [lindex \$args 0] { numatoms { return 2824 } name { return 1BL8.pdb } } }; return "" }
proc atomselect {args} { return ::sel }
proc ::sel {cmd args} { return "" }
proc measure {what sel} { return 20.607 }
namespace eval ::VMDPathFinder { variable state; array set state {} }
$(awk '/^proc ::VMDPathFinder::_mem_sig_field/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_run_dir_signature/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_discover_roots/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_run_structure_mismatch/,/^}/' "$TCL")
namespace eval ::VMDPathFinder {
    puts "IDS [join [lmapish [_mem_discover_roots "$TMP/hole_output_1BL8"]] ,]"
    puts "DIRS [join [dirs [_mem_discover_roots "$TMP/hole_output_1BL8"]] ,]"
    puts "IDS_PLAIN [join [lmapish [_mem_discover_roots "$TMP/plain"]] ,]"
    puts "SAME '[_run_structure_mismatch "$TMP/plain" 0]'"
    puts "OTHER '[_run_structure_mismatch "$TMP/other" 0]'"
}
TCLEOF
)
get() { printf '%s\n' "$OUT" | sed -n "s/^$1 //p" | head -1; }

[ "$(get IDS)" = "1,2,3" ] && ok "three runs under one folder become memories 1, 2, 3" || bad "discovered: $(get IDS)"
[ "$(get DIRS)" = "1BL8_20260910-120000-aaaaaa,1BL8_20260910-130000-bbbbbb,1BL8_20260910-140000-cccccc" ] \
  && ok "...in name (date) order, and a folder without frames is not one" || bad "dirs: $(get DIRS)"
[ "$(get IDS_PLAIN)" = "1" ] && ok "a run folder chosen directly is memory 1" || bad "plain: $(get IDS_PLAIN)"
[ "$(get SAME)" = "''" ] && ok "a run made on this structure loads without a word" || bad "same: $(get SAME)"
printf '%s' "$(get OTHER)" | grep -q "4649 atoms" && printf '%s' "$(get OTHER)" | grep -q "2824 atoms" \
  && ok "a run made on another structure is named, with both atom counts" || bad "other: $(get OTHER)"

echo "memory-load: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || { echo "--- raw ---"; printf '%s\n' "$OUT"; exit 1; }
