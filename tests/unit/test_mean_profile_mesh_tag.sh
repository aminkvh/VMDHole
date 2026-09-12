#!/bin/sh
# The Mean Profile cache filename must change when Settings > Surface mesher
# does, or a switch never rebuilds it - the OTHER mesher's tube keeps being
# served forever (neither plot_data_version nor _geomver moves on a mesher
# switch). build_and_show_mean_surface tags its files with [surface_mesh_tag 1]
# (union=1, matching the surface_mesh call it makes) for exactly this reason.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "mean-profile-mesh-tag"
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
[ -f "$TCL" ] || { echo "SKIP: no plugin source"; exit 0; }

OUT=$(tclsh <<TCLEOF 2>&1
namespace eval ::VMDPathFinder {
    variable state
    array set state {pore_method circular mesher csg csg_voxel 1.4 csg_voxel_fine 0.7}
    proc tool_path {name} { return "/bin/true" }
}
$(awk '/^proc ::VMDPathFinder::_csg_can_mesh/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_csg_voxel /,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::_csg_voxel_spec/,/^}/' "$TCL")
$(awk '/^proc ::VMDPathFinder::surface_mesh_tag/,/^}/' "$TCL")
namespace eval ::VMDPathFinder {
    puts "csg [surface_mesh_tag 1]"
    set state(mesher) sos
    puts "sos [surface_mesh_tag 1]"
    set state(mesher) csg
    set state(csg_voxel) 0.5
    puts "csg_finer [surface_mesh_tag 1]"
}
TCLEOF
)
csg_tag=$(printf '%s\n' "$OUT" | awk '$1=="csg"{print $2}')
sos_tag=$(printf '%s\n' "$OUT" | awk '$1=="sos"{print $2}')
finer_tag=$(printf '%s\n' "$OUT" | awk '$1=="csg_finer"{print $2}')
[ -n "$csg_tag" ] && ok "the csg mesher produces a non-empty tag ($csg_tag)" || bad "csg tag empty: $OUT"
[ "$sos_tag" = "" ] && ok "the sos mesher produces the untagged (empty) name" || bad "sos tag: '$sos_tag'"
[ "$csg_tag" != "$sos_tag" ] && ok "switching mesher changes the tag, so the cache filename changes" || bad "tags matched: '$csg_tag'"
[ -n "$finer_tag" ] && [ "$finer_tag" != "$csg_tag" ] && ok "changing the csg voxel size also changes the tag" || bad "finer tag: '$finer_tag' vs '$csg_tag'"
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
