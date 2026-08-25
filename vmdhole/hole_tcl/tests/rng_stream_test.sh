#!/bin/sh
# GATE ZERO: does the Tcl RNG draw the same stream as the Fortran one?
#
# Nothing else in this port can be validated until this passes. HOLE is Monte
# Carlo - its own reruns differ by 0.0053 A - so without an identical draw
# sequence, a radius disagreement is indistinguishable from noise and a real
# porting error cannot be told from a correct one.
#
# holcal.f calls DRAND, which is NOT in the Apache source release. Recovered
# from the binary instead: `nm ~/hole2/exe/hole` puts drand_ at 0x22180, and
# disassembling it shows it seeds once from cseed (forcing an odd value) then
# calls _gfortran_rand - gfortran's RAND intrinsic, i.e. Park-Miller.
#
# Tolerance is 1e-6, not exact equality: gfortran's RAND returns a REAL*4, so
# the reference values are the true sequence rounded to single precision. The
# underlying integer recurrence is what has to match, and at 1e-7 agreement it
# demonstrably does.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
echo "=============================================================="
echo "rng-stream: Tcl Park-Miller vs gfortran RAND"

command -v gfortran >/dev/null 2>&1 || { echo "  SKIP  no gfortran - cannot build the reference"; exit 0; }

# BOTH seeds. 12345 is a generic one; 1 is checked separately because it is
# the DEFAULT, and gfortran's generator does not treat it like any other seed -
# a port that matches at 12345 can still be wrong on essentially every real run.
SEEDS="12345 1"
fail=0
for SEED in $SEEDS; do
cat > "$TMP/ref.f" <<EOF
      PROGRAM RNGREF
      IMPLICIT NONE
      INTEGER I
      REAL R
      CALL SRAND( $SEED)
      DO 10 I = 1, 50
        R = RAND(0)
        WRITE(*,'(F20.16)') R
10    CONTINUE
      END
EOF
gfortran -std=legacy -o "$TMP/ref" "$TMP/ref.f" 2>/dev/null || {
    echo "  SKIP  reference would not build"; exit 0; }
"$TMP/ref" > "$TMP/fortran.txt" 2>/dev/null

cat > "$TMP/tcl.tcl" <<EOF
source "$DIR/../hole.tcl"
hole::rng::srand $SEED
for {set i 0} {\$i < 50} {incr i} { puts [format "%20.16f" [hole::rng::rand]] }
EOF
tclsh "$TMP/tcl.tcl" > "$TMP/tcl.txt" 2>&1 || { echo "  FAIL  Tcl side errored (seed $SEED)"; cat "$TMP/tcl.txt"; exit 1; }

echo "  --- seed $SEED"
paste "$TMP/fortran.txt" "$TMP/tcl.txt" | awk '
{ d = $1 - $2; if (d < 0) d = -d; if (d > max) max = d; n++;
  if (d > 1e-6) bad++ }
END {
  printf "  draws compared: %d\n", n
  printf "  max difference: %.3e\n", max
  if (n < 50)    { print "  FAIL  short stream"; exit 1 }
  if (bad > 0)   { printf "  FAIL  %d draw(s) over 1e-6\n", bad; exit 1 }
  print "  PASS  the Tcl stream matches gfortran RAND draw for draw"
}' || fail=1
done
exit $fail
