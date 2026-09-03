#!/bin/sh
# Contract test for sos_triangle_fast --ionflow-project, the C form of the Ion
# Flow tab's water pass (VMDHole's _ion_flow_scan). The Tcl loop it replaces is
# the reference: for each candidate point, offset from the frame COM, min-image
# per box dimension, project on the axis (z), perpendicular distance (R), and
# the signed distance to the nearest sphere SURFACE (d3). Points at R >= scan_r
# are dropped; frames and points keep input order.
#
# THE DEFECT CLASS THIS GUARDS: a helper that returns plausible numbers with
# the wrong convention - no min-imaging, unsigned d3, or R measured from the
# COM instead of the axis - would colour water occupancy quietly wrong, since
# d3 is only ever compared with a shell threshold downstream. So every number
# is checked against an independent awk re-derivation, not just "it ran".
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/native/sos_triangle_fast.c"
VOR="$ROOT/native/voronoi"
CC="${CC:-cc}"
pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
echo "ionflow-project: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no sos_triangle_fast.c at $SRC"; exit 0; }
command -v "$CC" >/dev/null 2>&1 || { echo "SKIP: no C compiler ($CC)"; exit 0; }
command -v awk >/dev/null 2>&1 || { echo "SKIP: no awk"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM
if ! "$CC" -O2 -o "$T/sos" "$SRC" "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" -lm -lpthread 2>"$T/cc.log"; then
    echo "SKIP: build failed:"; head -5 "$T/cc.log"; exit 0
fi
"$T/sos" --hole-features </dev/null 2>/dev/null | grep -q ionflowproject && ok "--hole-features advertises ionflowproject" || bad "ionflowproject missing from --hole-features"

# Two sphere sets, two frames; the axis is tilted so R is NOT a simple radial
# distance in the box frame, and box lengths make min-imaging bite (point 3 in
# frame 0 sits one box length away from the COM and must fold back).
cat > "$T/in.txt" <<'IN'
scan_r 6.0
S 0 3
0 0 0 2.0
0 0 4 1.5
0 0 8 3.0
F 0 0 0.5 0.25 -0.125 40 40 60 0.6 0 0.8 5
10 0.5 0.25 1.0
11 2.5 0.25 -0.125
12 0.5 20.25 -0.125
13 40.5 0.25 59.875
14 1.5 1.5 3.0
S 1 1
1 1 1 0.5
F 7 1 0 0 0 0 0 0 0 0 1 2
20 1 1 3
21 9 0 0
IN
"$T/sos" --ionflow-project "$T/in.txt" "$T/out.txt" 2>"$T/err.txt" || bad "exit status ($(head -1 "$T/err.txt"))"
[ -s "$T/out.txt" ] && ok "wrote output" || bad "no output"

# Independent re-derivation in awk (same conventions as the plugin's Tcl loop).
awk -v scan_r=6.0 '
function rnd(v){ return (v>=0) ? int(v+0.5) : -int(-v+0.5) }
function abs(v){ return v<0 ? -v : v }
/^S /{ id=$2; n=$3; ns[id]=n; for(k=0;k<n;k++){ getline; sx[id,k]=$1; sy[id,k]=$2; sz[id,k]=$3; sr[id,k]=$4 } ; next }
/^F /{ fr=$2; id=$3; cx=$4; cy=$5; cz=$6; Lx=$7; Ly=$8; Lz=$9; ux=$10; uy=$11; uz=$12; n=$13
       kept=0; out=""
       for(k=0;k<n;k++){ getline; idx=$1; rx=$2-cx; ry=$3-cy; rz=$4-cz
         if(Lx>0) rx-=Lx*rnd(rx/Lx); if(Ly>0) ry-=Ly*rnd(ry/Ly); if(Lz>0) rz-=Lz*rnd(rz/Lz)
         z=rx*ux+ry*uy+rz*uz; qx=rx-z*ux; qy=ry-z*uy; qz=rz-z*uz; R=sqrt(qx*qx+qy*qy+qz*qz)
         if(R>=scan_r) continue
         wx=cx+rx; wy=cy+ry; wz=cz+rz; best=1e30
         for(m=0;m<ns[id];m++){ dx=wx-sx[id,m]; dy=wy-sy[id,m]; dz=wz-sz[id,m]; d=sqrt(dx*dx+dy*dy+dz*dz)-sr[id,m]; if(d<best) best=d }
         out=out sprintf("%d %.9f %.9f %.9f\n", idx, z, R, best); kept++ }
       printf "F %d %d\n%s", fr, kept, out; next }' "$T/in.txt" > "$T/ref.txt"
# Compare at 1e-7 (the reference is awk double arithmetic; the binary prints 17 digits).
awk 'NR==FNR { ref[FNR]=$0; nref=FNR; next }
     { if (FNR>nref) { bad++; next }
       split(ref[FNR], a, " "); split($0, b, " ")
       if (a[1]!=b[1]) { bad++; next }
       if (a[1]=="F") { if (a[2]!=b[2] || a[3]!=b[3]) bad++; next }
       for (i=2;i<=4;i++) { d=a[i]-b[i]; if (d<0) d=-d; if (d>1e-7) bad++ }
       n++ }
     END { printf "%d %d %d\n", (bad+0), (n+0), (FNR==nref) }' "$T/ref.txt" "$T/out.txt" > "$T/cmp.txt"
read nbad npts samelen < "$T/cmp.txt"
[ "$nbad" -eq 0 ] && [ "$samelen" -eq 1 ] && ok "every kept point matches the awk re-derivation ($npts points, |d| <= 1e-7)" || bad "$nbad mismatching lines vs awk reference (samelen=$samelen)"
# Specific conventions, checked on their own so a regression names itself.
grep -q '^F 0 4$' "$T/out.txt" && ok "frame 0 drops the one point outside scan_r (5 in, 4 kept)" || bad "frame 0 kept count wrong: $(grep '^F 0' "$T/out.txt")"
# Point 13 is the COM displaced by exactly (Lx, 0, Lz): min-imaged it lands ON
# the COM, so z = R = 0 and d3 is the COM's own distance to sphere 0's surface
# (|COM| - 2 = -1.427...). Without min-imaging R would be ~40 and it would be dropped.
awk '$1==13 { if ($2 > -1e-9 && $2 < 1e-9 && $3 < 1e-9 && $4 < -1.42 && $4 > -1.43) f=1 } END { exit f?0:1 }' "$T/out.txt" && ok "point one box length away is min-imaged back onto the COM (z=R=0, d3=-1.427)" || bad "min-imaging: $(awk '$1==13' "$T/out.txt")"
awk '$1==10 { if ($4 < 0) f=1 } END { exit f?0:1 }' "$T/out.txt" && ok "a point inside a sphere has negative d3" || bad "signed d3: $(awk '$1==10' "$T/out.txt")"
grep -q '^F 7 1$' "$T/out.txt" && ok "second frame keeps input order and its own sphere set" || bad "frame 7: $(grep '^F 7' "$T/out.txt")"
# Malformed input must fail loudly, not emit a partial answer.
printf 'scan_r 5\nF 0 3 0 0 0 0 0 0 0 0 1 1\n1 0 0 0\n' > "$T/bad.txt"
if "$T/sos" --ionflow-project "$T/bad.txt" "$T/bad_out.txt" 2>/dev/null; then bad "undefined sphere set accepted"; else ok "undefined sphere set is refused (nonzero exit)"; fi
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
