#!/bin/sh
# Compare a hole.tcl profile against the real binary's, point for point.
#
# Usage: profile_vs_reference.sh PDB ?extra hole.tcl args?
#
# Runs ~/hole2/exe/hole and hole.tcl on the SAME cards and reports max and RMS
# radius difference plus the two bottlenecks. Prints numbers and does not
# judge them: what counts as agreement depends on which stage is being
# validated, and HOLE's own rerun noise is 0.0053 A.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$DIR/../../.." && pwd)
PDB=${1:?usage: profile_vs_reference.sh PDB ?extra args?}
shift
[ -x "$HOME/hole2/exe/hole" ] || { echo "SKIP: no reference binary"; exit 0; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cp "$PDB" "$TMP/in.pdb"
RADSRC="$ROOT/vmdhole/tests/fixtures/simple.rad"
[ -f "$RADSRC" ] || RADSRC="$ROOT/native/stock_build/hole2/rad/simple.rad"
cp "$RADSRC" "$TMP/s.rad"
# Short filenames on purpose: HOLE's card fields are fixed width and silently
# truncate a long path, which reads as "cannot open radius file".
printf 'coord in.pdb\nradius s.rad\ncvect 0 0 1\nsample 0.5\nendrad 8.0\n' > "$TMP/ctrl.inp"
( cd "$TMP" && "$HOME/hole2/exe/hole" < ctrl.inp > out.txt 2>&1 )

# The reference profile is the table after the cenxyz.cvec header; take the
# first two numeric columns (axial coordinate, radius).
awk '/cenxyz/{f=1;next} f && NF>=2 && ($1+0==$1) && ($2+0==$2) {printf "%.4f %.4f\n",$1,$2}' \
    "$TMP/out.txt" | sort -n > "$TMP/ref.dat"
if [ ! -s "$TMP/ref.dat" ]; then echo "FAIL: no reference profile parsed"; tail -5 "$TMP/out.txt"; exit 1; fi

tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb "$TMP/in.pdb" -rad "$TMP/s.rad" \
    -cvect "0 0 1" -sample 0.5 -endrad 8.0 -csv "$TMP/tcl.csv" "$@" > "$TMP/tcl.log" 2>&1 \
    || { echo "FAIL: hole.tcl errored"; cat "$TMP/tcl.log"; exit 1; }
awk -F, 'NR>1 {printf "%.4f %.4f\n",$1,$2}' "$TMP/tcl.csv" | sort -n > "$TMP/tcl.dat"

echo "=============================================================="
echo "profile-vs-reference: $(basename "$PDB")"
awk '
NR==FNR { rz[FNR]=$1; rr[FNR]=$2; nref=FNR; next }
{ tz[FNR]=$1; tr[FNR]=$2; nt=FNR }
END {
  # Nearest-axial-coordinate pairing: the two engines do not necessarily land
  # on identical slice positions, so index-for-index would compare unrelated
  # points.
  for (i=1;i<=nt;i++) {
    best=1e30; bj=0
    for (j=1;j<=nref;j++) { d=tz[i]-rz[j]; if(d<0)d=-d; if(d<best){best=d;bj=j} }
    if (bj && best<=0.30) { e=tr[i]-rr[bj]; if(e<0)e=-e
      if(e>max){max=e; maxz=tz[i]}; sum+=e*e; n++ }
  }
  for (j=1;j<=nref;j++) if (rmin=="" || rr[j]<rmin) { rmin=rr[j]; rminz=rz[j] }
  for (i=1;i<=nt;i++)   if (tmin=="" || tr[i]<tmin) { tmin=tr[i]; tminz=tz[i] }
  printf "  reference slices %d, tcl slices %d, paired %d\n", nref, nt, n
  printf "  bottleneck  reference %.4f A at %.3f\n", rmin, rminz
  printf "  bottleneck  hole.tcl  %.4f A at %.3f\n", tmin, tminz
  if (n>0) printf "  radius diff max %.4f A (at %.3f), RMS %.4f A\n", max, maxz, sqrt(sum/n)
  else     print  "  no slices paired within 0.30 A - the two grids do not overlap"
}' "$TMP/ref.dat" "$TMP/tcl.dat"
