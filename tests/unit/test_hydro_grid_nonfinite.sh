#!/bin/sh
# Regression: a non-finite coordinate in a --hydro-atoms or --hydro-sph
# sidecar must not cause UNDEFINED BEHAVIOUR in hydro_compute_sphere_h's
# atom/sphere grid.
#
# THE DEFECT THIS GUARDS (a sibling of test_sos_nonfinite_input.sh's - same
# NaN-cast hazard, different grid): hydro_compute_sphere_h bins atoms and
# looks up spheres with `(int)((v-origin)/cell)`, undefined for a NaN or
# out-of-range v. Unlike the .sos dot grid (a hash table - an out-of-range
# index just lands in a valid bucket), the atom-binning half of this grid is
# a DENSE array indexed directly by (ak*gy+aj)*gx+ai with no bounds check
# before the write, so the UB cast risked an out-of-bounds ghead[]/gnext[]
# write, not merely a wrong answer.
#
# Two things are asserted, and the second matters as much as the first: this
# codebase's whole justification is byte-identical output, so a fix that
# removed the UB by changing results for real (finite) coordinates would be
# worse than the bug.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/native/sos_triangle_fast.c"
VOR="$ROOT/native/voronoi"
CC="${CC:-cc}"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "hydro-grid-nonfinite: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no sos_triangle_fast.c at $SRC"; exit 0; }
command -v "$CC" >/dev/null 2>&1 || { echo "SKIP: no C compiler ($CC)"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# ---- inputs -------------------------------------------------------------
# A tiny valid .sos (type 4 = surface point) so hydro colouring has triangles
# to walk; the grid hazard fires while loading the sidecars, before any
# triangle is touched, so the surface itself just needs to be non-empty.
{ i=0; while [ $i -lt 30 ]; do
    printf '4 %d.5 %d.25 %d.75 0 0 0\n' "$i" "$((i%7))" "$((i%5))"; i=$((i+1)); done
} > "$T/valid.sos"

# --hydro-sph reads fixed-column PDB (pdb_col: x=cols31-38, y=39-46, z=47-54,
# r=61-66, 1-indexed) - build with explicit substr placement, not printf
# field widths, so the NaN token lands in exactly the x column and nowhere
# else drifts.
sph_line() {
  awk -v s="$1" -v x="$2" -v y="$3" -v z="$4" -v r="$5" '
    BEGIN {
      line = sprintf("ATOM  %5d  QC1 SPH S%4d", s, s)
      while (length(line) < 66) line = line " "
      line = substr(line,1,30) sprintf("%8s",x) sprintf("%8s",y) sprintf("%8s",z) \
             substr(line,55,6) sprintf("%6s",r)
      print line
    }'
}
{ sph_line 1 nan  1.000 1.000 3.00
  sph_line 2 1.000 2.000 3.000 3.00
  sph_line 3 -3.000 -2.000 -1.000 2.50
  sph_line 4 4.000 0.500 2.000 2.50
} > "$T/nonfinite.sph"
{ sph_line 1 1.000 1.000 1.000 3.00
  sph_line 2 1.000 2.000 3.000 3.00
  sph_line 3 -3.000 -2.000 -1.000 2.50
  sph_line 4 4.000 0.500 2.000 2.50
} > "$T/valid.sph"

# --hydro-atoms reads plain "x y z hkd hww" text (strtod, not pdb_col) - a
# bare "nan" token parses the same way "1.0" does.
{ echo "nan 1.0 1.0 -0.5 -0.5"
  echo "1.0 1.0 1.0 -0.5 -0.5"
  echo "2.0 1.5 0.5 0.3 0.3"
  echo "-1.0 -1.0 -1.0 0.8 0.8"
  echo "3.0 2.0 1.0 -1.2 -1.2"
} > "$T/nonfinite_atoms.dat"
{ echo "1.0 1.0 1.0 -0.5 -0.5"
  echo "2.0 1.5 0.5 0.3 0.3"
  echo "-1.0 -1.0 -1.0 0.8 0.8"
  echo "3.0 2.0 1.0 -1.2 -1.2"
} > "$T/valid_atoms.dat"

# ---- build ----------------------------------------------------------------
if ! "$CC" -O1 -o "$T/sos" "$SRC" "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" \
        -lm -lpthread 2>"$T/build.log"; then
    echo "SKIP: could not build sos_triangle_fast"; sed 's/^/    /' "$T/build.log" | head -5; exit 0
fi

SANMODE="crash-only (no sanitizer)"
if "$CC" -O1 -fsanitize=undefined -fno-sanitize-recover=undefined -o "$T/sos_ub" "$SRC" \
        "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" -lm -lpthread 2>/dev/null; then
    SANMODE="UBSan"
fi
echo "  mode: $SANMODE"
BIN="$T/sos"; [ "$SANMODE" = "UBSan" ] && BIN="$T/sos_ub"

# ---- 1. a non-finite ATOM coordinate (--hydro-atoms) must not trip UB -----
UBSAN_OPTIONS=print_stacktrace=0 "$BIN" -s --hydro-atoms "$T/nonfinite_atoms.dat" \
    --hydro-sph "$T/valid.sph" --hydro-scheme kd \
    < "$T/valid.sos" > "$T/nf_atoms.out" 2> "$T/nf_atoms.err"
rc=$?
if [ "$rc" -ge 128 ]; then
    bad "non-finite hydro atom killed the process with signal $((rc-128))"
elif grep -q 'runtime error:' "$T/nf_atoms.err"; then
    bad "non-finite hydro atom triggers undefined behaviour"
    grep 'runtime error:' "$T/nf_atoms.err" | head -2 | sed 's/^/          /'
else
    ok "non-finite hydro atom handled without UB or a signal (exit $rc)"
fi

# ---- 2. a non-finite SPHERE coordinate (--hydro-sph) must not trip UB -----
UBSAN_OPTIONS=print_stacktrace=0 "$BIN" -s --hydro-atoms "$T/valid_atoms.dat" \
    --hydro-sph "$T/nonfinite.sph" --hydro-scheme kd \
    < "$T/valid.sos" > "$T/nf_sph.out" 2> "$T/nf_sph.err"
rc=$?
if [ "$rc" -ge 128 ]; then
    bad "non-finite hydro sphere killed the process with signal $((rc-128))"
elif grep -q 'runtime error:' "$T/nf_sph.err"; then
    bad "non-finite hydro sphere triggers undefined behaviour"
    grep 'runtime error:' "$T/nf_sph.err" | head -2 | sed 's/^/          /'
else
    ok "non-finite hydro sphere handled without UB or a signal (exit $rc)"
fi

# ---- 3. valid sidecars must be BYTE-IDENTICAL and non-empty across builds -
# The point is that a future edit to the grid cast cannot silently move a
# real, all-finite colouring - not just that this run happens to not crash.
"$T/sos" -s --hydro-atoms "$T/valid_atoms.dat" --hydro-sph "$T/valid.sph" \
    --hydro-scheme kd < "$T/valid.sos" > "$T/v1.out" 2>/dev/null
"$T/sos_ub" -s --hydro-atoms "$T/valid_atoms.dat" --hydro-sph "$T/valid.sph" \
    --hydro-scheme kd < "$T/valid.sos" > "$T/v2.out" 2>/dev/null
if [ "$SANMODE" != "UBSan" ]; then
    "$T/sos" -s --hydro-atoms "$T/valid_atoms.dat" --hydro-sph "$T/valid.sph" \
        --hydro-scheme kd < "$T/valid.sos" > "$T/v2.out" 2>/dev/null
fi
if cmp -s "$T/v1.out" "$T/v2.out" && [ -s "$T/v1.out" ]; then
    ok "valid hydro sidecars produce stable, non-empty output ($(grep -c 'draw ' "$T/v1.out") draw records)"
else
    bad "valid hydro sidecars output is empty or the plain/UBSan builds disagree"
fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
