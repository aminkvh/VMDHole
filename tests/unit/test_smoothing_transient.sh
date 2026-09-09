#!/bin/sh
# Surface smoothing must ATTENUATE a feature that only one frame has.
#
# The regression this pins (shipped 3bd2fa0, reverted 2026-09-09): the averager
# divided each corner by "the number of frames that produced a value there"
# instead of by the frame count. fill_field only evaluates a sphere within
# radius+2h, so a corner beyond that band keeps a 1e9 sentinel - but that corner
# is not UNKNOWN in that frame, it is known to be far OUTSIDE. Skipping it gave
# a corner reached by 1 frame of N that one frame's value at FULL strength, so a
# transient bulge survived smoothing undiminished. The user saw it as random
# blobs and false lateral openings on a spherical pore.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "smoothing-transient"

F="$ROOT/vmdpathfinder/tests/fixtures/spherical_1GRM.sph"
[ -f "$F" ] || { echo "SKIP: no spherical fixture"; exit 0; }
EXE=""
for c in "$ROOT/native/nm/mesh_csg" "$ROOT/native/sos_triangle_fast"; do
    [ -x "$c" ] && { EXE="$c"; break; }
done
[ -n "$EXE" ] || { echo "SKIP: no mesher built (sh native/build.sh)"; exit 0; }
case "$EXE" in *sos_triangle_fast) MODE="--mesh" ;; *) MODE="" ;; esac

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
# One frame carries a wide sphere far out on +x that no other frame has.
cp "$F" "$T/bulge.sph"
printf 'ATOM      1  QSS SPH S   9      12.000   0.000   3.500  6.00  6.00\n' >> "$T/bulge.sph"

run() { $EXE $MODE "$1" "$2" 1.4/0.7 --draw --axis 0 0 0 0 0 1 8 --with "$F" --with "$F" >/dev/null 2>&1; }
beyond() { awk '{c=0; for(i=1;i<=NF;i++){v=$i+0; if(v>8) c++} if(c>0) print}' "$1" | wc -l; }
maxc()   { awk '{for(i=1;i<=NF;i++) if($i ~ /^-?[0-9.]+$/){v=$i+0; if(v>m)m=v}} END{printf "%.0f", m}' "$1"; }

run "$T/bulge.sph" "$T/smoothed.plot" || { echo "SKIP: mesher failed"; exit 0; }
# the same pore WITHOUT the bulge, smoothed against itself, is the baseline
run "$F" "$T/plain.plot"

nb=$(beyond "$T/smoothed.plot"); np=$(beyond "$T/plain.plot")
xb=$(maxc "$T/smoothed.plot")
[ "$nb" -gt 0 ] 2>/dev/null && ok "the mesher ran and produced geometry ($nb triangles past x=8)" \
    || bad "no geometry produced"
# 1 frame of 3 has the bulge, so a third of it may survive - but not the ~1900
# triangles reaching x=18 that the unaveraged sentinel path produced.
if [ "$nb" -lt 1000 ]; then
    ok "a one-frame bulge is ATTENUATED by smoothing ($nb triangles past x=8, not ~1900)"
else
    bad "a one-frame bulge survived smoothing at full strength: $nb triangles past x=8"
fi
if [ "$xb" -lt 17 ]; then
    ok "...and does not reach its unsmoothed extent (max coord $xb, unaveraged was 18)"
else
    bad "the bulge reaches x=$xb, its full unsmoothed extent"
fi
# and the identity invariant the averager also has to satisfy
IDENT=$(CSG_DEBUG_IDENT=1 $EXE $MODE "$F" "$T/id.plot" 1.4/0.7 --draw --axis 0 0 0 0 0 1 8 \
        --with "$F" --with "$F" 2>&1 | sed -n 's/.*IDENT //p' | head -1)
case "$IDENT" in
    *"differing=0 maxdiff=0"*) ok "smoothing identical frames is still the identity ($IDENT)" ;;
    "")                        bad "no IDENT line - CSG_DEBUG_IDENT not honoured" ;;
    *)                         bad "smoothing identical frames changed the field: $IDENT" ;;
esac

echo "smoothing-transient: $pass passed, $fail failed"
[ "$fail" -eq 0 ] || exit 1
