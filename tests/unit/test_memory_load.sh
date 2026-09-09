#!/bin/sh
# Multi-memory load: discovery of the slots saved under one root, and the gate
# that refuses to show together analyses that were not computed the same way.
#
# The gate matters because load is the ONE path that can build the mixed state
# the shared settings otherwise make unreachable: memory 1 saved at sample 0.25
# next to memory 2 saved at 1.0 looks like a structural difference and is not.
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
mksig() { # dir sample method search
  mkdir -p "$1/frame_0"
  printf 'sel|protein|cpoint|1 2 3|sample|%s|endrad|15.0|method|%s|search|%s|conneng|hole|plugin|1.0' \
    "$2" "$3" "$4" > "$1/frame_0/vmdpathfinder_signature.dat"
}
# a coherent root: three memories, same sampling / method / engine
mksig "$TMP/good"        0.25 circular mc
mksig "$TMP/good/mem_2"  0.25 circular mc
mksig "$TMP/good/mem_3"  0.25 circular mc
# an incoherent root: memory 2 was run at a different sampling
mksig "$TMP/bad"         0.25 circular mc
mksig "$TMP/bad/mem_2"   1.0  circular mc
# a second incoherent root: different search engine
mksig "$TMP/bad2"        0.25 circular mc
mksig "$TMP/bad2/mem_2"  0.25 circular nm
# a folder with no mem_* at all - the ordinary single-run case
mksig "$TMP/plain"       0.25 circular mc
# a mem_* dir with no frames must not be offered as a memory
mkdir -p "$TMP/good/mem_9"
# memory 1 never run: the root holds ONLY mem_2
mksig "$TMP/only2/mem_2" 0.25 circular mc
mkdir -p "$TMP/only2"

# lmapish: Tcl 8.5 has no lmap, and the plugin targets 8.5 only.
OUT=$(tclsh <<TCLEOF 2>&1
proc lmapish {pairs} { set o {}; foreach p \$pairs { lappend o [lindex \$p 0] }; return \$o }
namespace eval ::VMDPathFinder { variable state; array set state {} }
$(awk '/^proc ::VMDPathFinder::_mem_sig_field/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_sig_shared_fields/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_run_dir_signature/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_discover_roots/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_load_incoherent/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_sync_shared_labels/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_mem_adopt_shared_from_sig/,/^}/' "$TCL")
namespace eval ::VMDPathFinder {
    puts "IDS_GOOD [join [lmapish [_mem_discover_roots "$TMP/good"]] ,]"
    puts "IDS_PLAIN [join [lmapish [_mem_discover_roots "$TMP/plain"]] ,]"
    puts "IDS_ONLY2 [join [lmapish [_mem_discover_roots "$TMP/only2"]] ,]"
    # the Load button routes to the memory loader whenever a slot other than 1
    # is present, even when that is the ONLY slot
    set _p [_mem_discover_roots "$TMP/only2"]
    set _need [expr {[llength \$_p] > 1}]
    foreach _q \$_p { if {[lindex \$_q 0] != 1} { set _need 1 } }
    puts "ROUTE_ONLY2 \$_need"
    set _p [_mem_discover_roots "$TMP/plain"]
    set _need [expr {[llength \$_p] > 1}]
    foreach _q \$_p { if {[lindex \$_q 0] != 1} { set _need 1 } }
    puts "ROUTE_PLAIN \$_need"
    puts "GOOD_BAD '[_mem_load_incoherent [_mem_discover_roots "$TMP/good"]]'"
    set m [_mem_load_incoherent [_mem_discover_roots "$TMP/bad"]]
    puts "BAD_SAMPLE [expr {\$m ne "" && [string match "*sample*" \$m]}]"
    set m2 [_mem_load_incoherent [_mem_discover_roots "$TMP/bad2"]]
    puts "BAD_SEARCH [expr {\$m2 ne "" && [string match "*search*" \$m2]}]"
    puts "MSG \$m"
    _mem_adopt_shared_from_sig [_mem_run_dir_signature "$TMP/bad/mem_2"]
    puts "ADOPT \$state(sample) \$state(pore_method) \$state(search_engine) \$state(pore_method_disp)"
}
TCLEOF
)
get() { printf '%s\n' "$OUT" | sed -n "s/^$1 //p" | head -1; }

[ "$(get IDS_GOOD)" = "1,2,3" ] && ok "all three saved memories are discovered under one root" || bad "discovered: $(get IDS_GOOD)"
printf '%s' "$(get IDS_GOOD)" | grep -q "9" && bad "an empty mem_9 was offered as a memory" || ok "a mem_* folder with no frames is not offered"
[ "$(get IDS_PLAIN)" = "1" ] && ok "an ordinary single-run folder is just memory 1" || bad "plain: $(get IDS_PLAIN)"
[ "$(get GOOD_BAD)" = "''" ] && ok "a coherent set loads" || bad "coherent set refused: $(get GOOD_BAD)"
[ "$(get BAD_SAMPLE)" = "1" ] && ok "different SAMPLING between memories is refused" || bad "sampling mismatch not caught"
[ "$(get BAD_SEARCH)" = "1" ] && ok "different SEARCH ENGINE between memories is refused" || bad "engine mismatch not caught"
printf '%s' "$(get MSG)" | grep -q "memory 1" && printf '%s' "$(get MSG)" | grep -q "memory 2" \
  && ok "the refusal names both memories" || bad "refusal message unhelpful: $(get MSG)"
[ "$(get IDS_ONLY2)" = "2" ] && ok "a root where memory 1 was never run still reports slot 2" || bad "only2: $(get IDS_ONLY2)"
[ "$(get ROUTE_ONLY2)" = "1" ] && ok "...and the Load button routes it to the memory loader, so it does not come back as memory 1" || bad "only2 took the plain load path"
[ "$(get ROUTE_PLAIN)" = "0" ] && ok "an ordinary single-run folder still takes the plain load path" || bad "plain folder diverted to the memory loader"
[ "$(get ADOPT)" = "1.0 circular mc Spherical" ] && ok "a load adopts the settings the run was computed with, label included" || bad "adopt: $(get ADOPT)"

echo "memory-load: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || { echo "--- raw ---"; printf '%s\n' "$OUT"; exit 1; }
