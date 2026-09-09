#!/bin/sh
# A capsule surface built through sph_process + sos_triangle must cover the
# whole pore, not a fraction of it.
#
# THE DEFECT THIS GUARDS:
#   sph_process's -colour mode also emits "endrad points" at colour -1, which
#   sos_triangle uses to cut sharp ends on a SPHERICAL pore. On capsule records
#   those markers clip nearly the entire surface away. The plugin passed
#   -colour unconditionally, so capsule surfaces rendered as a stub: measured on
#   one frame, 1596 triangles spanning 22.8 A where the marching-cubes mesher
#   drew 148.66 A from the same spheres. Reported as "sos in capsule does not
#   render the full surface but the marching does".
#
# This drives the REAL binaries on a committed capsule .sph, both ways, so it
# fails if either the flag's effect changes or the plugin starts passing it
# again.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SPH="$ROOT/vmdpathfinder/tests/fixtures/capsule_1GRM.sph"
TCL="$ROOT/vmdpathfinder/vmdpathfinder.tcl"
SP="${SPH_PROCESS:-$HOME/hole2/exe/sph_process}"
ST="${SOS_TRIANGLE:-$HOME/hole2/exe/sos_triangle}"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "capsule-sos-full-surface"
[ -f "$SPH" ] || { echo "SKIP: no capsule .sph fixture"; exit 0; }
[ -x "$SP" ] && [ -x "$ST" ] || { echo "SKIP: sph_process/sos_triangle not available"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

span() {   # $1 = .vmd_plot -> axial span of its triangles, or empty
    awk '/trinorm/{gsub(/[{}]/," "); c=0
           for(i=1;i<=NF;i++){ if($i ~ /^-?[0-9.]+$/){ c++
               if(c%3==0){ z=$i+0; if(m++==0){lo=hi=z} if(z<lo)lo=z; if(z>hi)hi=z } }
               if(c>=9) break } }
         END{ if(m) printf "%.2f", hi-lo }' "$1"
}
build() {  # $1 = extra sph_process flags, $2 = tag
    "$SP" -sos -dotden 15 $1 "$SPH" "$T/$2.sos" >/dev/null 2>&1 || return 1
    "$ST" -s < "$T/$2.sos" > "$T/$2.plot" 2>/dev/null || return 1
    return 0
}

# the .sph's own extent is the yardstick
SPHSPAN=$(awk 'substr($0,1,4)=="ATOM"{z=substr($0,47,8)+0; if(n++==0){lo=hi=z} if(z<lo)lo=z; if(z>hi)hi=z}
               END{printf "%.2f", hi-lo}' "$SPH")

if build "" plain && [ -n "$(span "$T/plain.plot")" ]; then
    PLAIN=$(span "$T/plain.plot")
    # The mesh envelopes the spheres, so its extent must be AT LEAST theirs.
    # Not a fraction of it: a 0.8 factor passed trivially here (58.8 vs 20.1)
    # and would have passed on the truncated build too if the pore were shorter.
    # The truncated build measures 4.31 A against a 20.11 A sphere extent, so
    # this comparison is what actually catches the regression.
    if awk -v a="$PLAIN" -v b="$SPHSPAN" 'BEGIN{exit !(a >= b)}'; then
        ok "capsule surface covers the pore without -colour ($PLAIN of $SPHSPAN A)"
    else
        bad "capsule surface is short without -colour ($PLAIN of $SPHSPAN A)"
    fi
else
    bad "capsule surface could not be built without -colour"
fi

# The hazard itself: -colour must be shown to truncate, or this test proves
# nothing about why the plugin omits it.
if build "-colour" col && [ -n "$(span "$T/col.plot")" ]; then
    COL=$(span "$T/col.plot")
    if awk -v a="$COL" -v b="$PLAIN" 'BEGIN{exit !(a < 0.8*b)}'; then
        ok "-colour still truncates it, which is why the plugin omits it ($COL vs $PLAIN A)"
    else
        ok "-colour no longer truncates ($COL vs $PLAIN A) - the workaround may be droppable"
    fi
fi

# ...and the plugin must actually omit it, in BOTH surface builders.
n=$(grep -c 'set color 0' "$TCL" 2>/dev/null || echo 0)
if [ "$n" -ge 2 ]; then
    ok "both surface builders drop -colour under capsule"
else
    bad "expected 2 capsule 'set color 0' guards in the surface builders, found $n"
fi

echo "capsule-sos-full-surface: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
