#!/bin/sh
# The Connolly cloud is thinned before triangulation only where it is cheap
# and invisible: under the marching-cubes mesher and in the playback draft.
# sph_process + sos_triangle draws the exact sphere union, and the thinned
# cloud shows as beads, so the settle pass keeps the full cloud.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "conn-render-target"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

OUT=$(tclsh <<TCLEOF 2>&1
namespace eval ::VMDPathFinder {
    variable state
    array set state {pore_method connolly mesher csg}
    proc tool_path {name} { return "/bin/true" }
}
$(awk '/^proc ::VMDPathFinder::_csg_can_mesh/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_conn_render_target_npts/,/^}/' "$TCL")
namespace eval ::VMDPathFinder {
    puts "csg [_conn_render_target_npts]"
    puts "csg_draft [_conn_render_target_npts 1]"
    set state(mesher) sos
    puts "sos [_conn_render_target_npts]"
    puts "sos_draft [_conn_render_target_npts 1]"
}
TCLEOF
)
case "$OUT" in *"csg 14000"*)       ok "marching-cubes mesher thins to 14000" ;;      *) bad "csg: $OUT" ;; esac
case "$OUT" in *"csg_draft 14000"*) ok "draft under the mesher thins to 14000" ;;    *) bad "csg draft: $OUT" ;; esac
case "$OUT" in *"sos 120000"*)      ok "sos settle pass keeps the full cloud" ;;     *) bad "sos: $OUT" ;; esac
case "$OUT" in *"sos_draft 14000"*) ok "sos draft still thins to 14000" ;;           *) bad "sos draft: $OUT" ;; esac
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
