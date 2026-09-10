#!/bin/sh
# A capsule surface built through sph_process + sos_triangle must cover the
# whole pore, not a fraction of it.
#
# sph_process -colour files a capsule's third radius band under the end-cap
# header, which sos_triangle clips: the surface came out as a stub. The plugin
# repairs that header after sph_process. This drives the real binaries on a
# committed capsule .sph.
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

# The hazard: -colour files the third radius band under the end-cap header
# {1 -1 -1 -1}, which sos_triangle clips. The plugin keeps -colour and repairs
# that header after sph_process (_capsule_sos_fix / _capsule_sos_fix_cmd), so
# the coloured surface must span the pore like the plain one.
if build "-colour" col && [ -n "$(span "$T/col.plot")" ]; then
    COL=$(span "$T/col.plot")
    awk 'NF==7 && $1+0==1 && $2+0==-1 && $3+0==-1 && $4+0==-1 {printf "%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n",1,2,-55,18,0,0,0; next} {print}' \
        "$T/col.sos" > "$T/fix.sos"
    "$ST" -s < "$T/fix.sos" > "$T/fix.plot" 2>/dev/null
    FIX=$(span "$T/fix.plot")
    NCOL=$(grep -o "color [a-z0-9]*" "$T/fix.plot" | sort -u | wc -l)
    if awk -v a="$FIX" -v b="$SPHSPAN" 'BEGIN{exit !(a >= b)}'; then
        ok "with the header repaired, the coloured surface covers the pore ($FIX of $SPHSPAN A; unrepaired $COL A)"
    else
        bad "the repaired coloured surface is still short ($FIX of $SPHSPAN A)"
    fi
    [ "$NCOL" -ge 3 ] && ok "...in three radius bands" || bad "expected three colour bands, found $NCOL"
else
    bad "capsule surface could not be built with -colour"
fi

# ...and the plugin repairs it in every builder that runs sph_process.
n=$(grep -c '_capsule_sos_fix' "$TCL" 2>/dev/null || echo 0)
if [ "$n" -ge 4 ]; then
    ok "the plugin keeps -colour and repairs the header in its builders ($n sites)"
else
    bad "expected the header repair in the surface builders, found $n mention(s)"
fi

echo "capsule-sos-full-surface: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
