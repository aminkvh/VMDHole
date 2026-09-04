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
# ---- the streamed form ------------------------------------------------------
# The plugin's fast path does not send per-frame filtered points at all: it
# sends every water coordinate as a raw little-endian float32 stream and lets
# the binary apply the cylinder and the axial window. That path must agree
# EXACTLY with the inline one it replaces - it is the same arithmetic on the
# same numbers, and a silent divergence would colour water occupancy wrongly
# with nothing to notice it by. Every coordinate below is exact in float32, so
# "exactly" means bit-for-bit, not "close".
cat > "$T/mkbin.c" <<'C'
/* Reads "x y z" triples and writes them as the stream layout the binary
   expects: per frame, nw x values, then nw y, then nw z, little-endian
   float32. Frames are separated by a blank line. */
#include <stdio.h>
#include <string.h>
static void put_le(FILE *o, double v) {
    float f = (float)v; unsigned int u; unsigned char b[4];
    memcpy(&u, &f, 4);
    b[0]=(unsigned char)(u&0xff); b[1]=(unsigned char)((u>>8)&0xff);
    b[2]=(unsigned char)((u>>16)&0xff); b[3]=(unsigned char)((u>>24)&0xff);
    fwrite(b,1,4,o);
}
int main(void) {
    char line[256]; double x[64],y[64],z[64]; int n=0,i;
    FILE *o = stdout;
    while (fgets(line,sizeof line,stdin)) {
        if (sscanf(line,"%lf %lf %lf",&x[n],&y[n],&z[n]) == 3) { n++; continue; }
        if (n) { for(i=0;i<n;i++) put_le(o,x[i]);
                 for(i=0;i<n;i++) put_le(o,y[i]);
                 for(i=0;i<n;i++) put_le(o,z[i]); n=0; }
    }
    if (n) { for(i=0;i<n;i++) put_le(o,x[i]);
             for(i=0;i<n;i++) put_le(o,y[i]);
             for(i=0;i<n;i++) put_le(o,z[i]); }
    return 0;
}
C
if ! "$CC" -O1 -o "$T/mkbin" "$T/mkbin.c" 2>"$T/cc2.log"; then
    echo "SKIP: helper build failed:"; head -3 "$T/cc2.log"
else
# Two frames, five points each, the same sphere set - written once as inline
# points and once as a stream, both with the same axial window.
POINTS_A='0.5 0.25 1.0
2.5 0.25 -0.125
0.5 20.25 -0.125
40.5 0.25 59.875
1.5 1.5 3.0'
POINTS_B='1.0 1.0 3.0
9.0 0.0 0.0
0.5 0.25 1.0
-3.5 2.0 0.5
2.5 0.25 -0.125'
{
  echo "scan_r 6.0"
  echo "zwin -20 20"
  echo "S 0 3"
  echo "0 0 0 2.0"; echo "0 0 4 1.5"; echo "0 0 8 3.0"
  echo "F 0 0 0.5 0.25 -0.125 40 40 60 0.6 0 0.8 5"
  echo "$POINTS_A" | awk '{printf "%d %s %s %s\n", NR+9, $1, $2, $3}'
  echo "F 7 0 0.5 0.25 -0.125 40 40 60 0.6 0 0.8 5"
  echo "$POINTS_B" | awk '{printf "%d %s %s %s\n", NR+9, $1, $2, $3}'
} > "$T/inline.txt"
printf '%s
' "$POINTS_A" "" "$POINTS_B" | "$T/mkbin" > "$T/stream.bin"
{
  echo "scan_r 6.0"
  echo "zwin -20 20"
  echo "coords 5"
  echo "$T/stream.bin"
  echo "group 1"
  echo "index 5"
  seq 10 14
  echo "S 0 3"
  echo "0 0 0 2.0"; echo "0 0 4 1.5"; echo "0 0 8 3.0"
  echo "F 0 0 0.5 0.25 -0.125 40 40 60 0.6 0 0.8 -1"
  echo "F 7 0 0.5 0.25 -0.125 40 40 60 0.6 0 0.8 -1"
} > "$T/stream.txt"
"$T/sos" --ionflow-project "$T/inline.txt" "$T/inline_out.txt" 2>/dev/null || bad "inline job failed"
"$T/sos" --ionflow-project "$T/stream.txt" "$T/stream_out.txt" 2>/dev/null || bad "streamed job failed"
grep -q '^T ' "$T/stream_out.txt" && ok "streamed job emits grouped T records" || bad "no T records: $(head -2 "$T/stream_out.txt")"
# Flatten both to sorted "idx frame z R d3" tuples and compare byte for byte.
awk '/^F /{fr=$2; next} {print $1, fr, $2, $3, $4}' "$T/inline_out.txt" | sort > "$T/flat_inline.txt"
awk '/^T /{idx=$2; n=$3; for(i=1;i<=n;i++){f[i]=$(3+i); z[i]=$(3+n+i); r[i]=$(3+2*n+i); d[i]=$(3+3*n+i)}
           for(i=1;i<=n;i++) print idx, f[i], z[i], r[i], d[i]}' "$T/stream_out.txt" | sort > "$T/flat_stream.txt"
if cmp -s "$T/flat_inline.txt" "$T/flat_stream.txt"; then
    ok "streamed output is bit-identical to the inline output ($(wc -l < "$T/flat_stream.txt") samples over 2 frames)"
else
    bad "streamed and inline disagree:"; diff "$T/flat_inline.txt" "$T/flat_stream.txt" | head -6
fi
[ -s "$T/flat_stream.txt" ] && ok "the shared window kept some samples (an empty match would prove nothing)" || bad "both outputs are empty"
# The axial window must actually bite: point 13 min-images onto the COM (z=0,
# inside the window) while a point beyond it must be dropped.
awk '$1==13' "$T/flat_stream.txt" | grep -q . && ok "min-imaged point survives the streamed path too" || bad "point 13 missing from the streamed output"
sed 's/^zwin -20 20$/zwin 0.5 20/' "$T/stream.txt" > "$T/stream_win.txt"
"$T/sos" --ionflow-project "$T/stream_win.txt" "$T/win_out.txt" 2>/dev/null || bad "windowed job failed"
if awk '/^T /{idx=$2; n=$3; for(i=1;i<=n;i++) if ($(3+n+i) <= 0.5) bad=1} END{exit bad?1:0}' "$T/win_out.txt"; then
    ok "a tightened axial window drops every sample below it"
else
    bad "a sample below the axial window survived"
fi
# Frames whose points come from the stream need the stream to be there.
sed 's/^group 1$//' "$T/stream.txt" > "$T/nogroup.txt"
"$T/sos" --ionflow-project "$T/nogroup.txt" "$T/nogroup_out.txt" 2>/dev/null \
  && grep -q '^F ' "$T/nogroup_out.txt" && ok "a streamed job without 'group' falls back to per-frame output" \
  || bad "streamed job without 'group' did not produce per-frame output"
grep -v '^index 5$' "$T/stream.txt" | grep -vE '^1[0-4]$' > "$T/noindex.txt"
if "$T/sos" --ionflow-project "$T/noindex.txt" "$T/noindex_out.txt" 2>/dev/null; then
    bad "a streamed job with no index list was accepted"
else
    ok "a streamed job with no index list is refused"
fi
fi

# Malformed input must fail loudly, not emit a partial answer.
printf 'scan_r 5\nF 0 3 0 0 0 0 0 0 0 0 1 1\n1 0 0 0\n' > "$T/bad.txt"
if "$T/sos" --ionflow-project "$T/bad.txt" "$T/bad_out.txt" 2>/dev/null; then bad "undefined sphere set accepted"; else ok "undefined sphere set is refused (nonzero exit)"; fi
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
