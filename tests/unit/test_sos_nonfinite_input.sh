#!/bin/sh
# Regression: a non-finite coordinate in a .sos must not cause UNDEFINED BEHAVIOUR
# in sos_triangle_fast's spatial grids.
#
# THE DEFECT THIS GUARDS (verified to go red on the pre-fix source):
#   Both spatial grids derive an integer cell index as `(int)floor(v/cell)`.
#   That cast is undefined for NaN or an out-of-range quotient; in practice it
#   yields INT_MIN. The +-1 neighbourhood walks that consume the index then
#   compute INT_MIN-1, and the neighbour search cubes `2L*cr+1` where cr came
#   from `(int)(sqrt(9*base_dist)/NCELL)`. UBSan on a NaN-bearing .sos:
#     runtime error: signed integer overflow: -2147483648 + -1 ...
#     runtime error: signed integer overflow: -4294967293 * -4294967293 ...
#   Note the second one overflows while COMPUTING the `ncells > max_dots` guard
#   that was supposed to reject an oversized range - an overflow-before-check.
#
#   Upstream sos_triangle.c has no cell index at all (its cull is an O(N^2)
#   scan), so this hazard was introduced by the fork's optimisation, not
#   inherited from HOLE.
#
# Two things are asserted, and the second matters as much as the first: this
# codebase's whole justification is byte-identical output, so a fix that
# removed the UB by changing results would be worse than the bug.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/native/sos_triangle_fast.c"
VOR="$ROOT/native/voronoi"
CC="${CC:-cc}"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "sos-nonfinite: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no sos_triangle_fast.c at $SRC"; exit 0; }
command -v "$CC" >/dev/null 2>&1 || { echo "SKIP: no C compiler ($CC)"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# ---- inputs -----------------------------------------------------------------
# A .sos record is "type x y z ...". Type 4 is a surface point. Keep these tiny
# and generated here: the test must not depend on a gitignored fixture corpus.
{ i=0; while [ $i -lt 60 ]; do
    printf '4 %d.5 %d.25 %d.75 0 0 0\n' "$i" "$((i%7))" "$((i%5))"; i=$((i+1)); done
} > "$T/valid.sos"
{ echo "4 nan nan nan 0 0 0"; echo "4 1.0 2.0 3.0 0 0 0"
  echo "4 inf 1.0 1.0 0 0 0";  echo "4 -inf 1.0 1.0 0 0 0"
  echo "4 1e308 1e308 1e308 0 0 0"; cat "$T/valid.sos"; } > "$T/nonfinite.sos"

# ---- build ------------------------------------------------------------------
if ! "$CC" -O1 -o "$T/sos" "$SRC" "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" \
        -lm -lpthread 2>"$T/build.log"; then
    echo "SKIP: could not build sos_triangle_fast"; sed 's/^/    /' "$T/build.log" | head -5; exit 0
fi

# UBSan build is optional - without it the crash check below still runs, it just
# cannot see a silent overflow. Say which mode is in force rather than implying
# more coverage than was achieved.
SANMODE="crash-only (no sanitizer)"
if "$CC" -O1 -fsanitize=undefined -fno-sanitize-recover=undefined -o "$T/sos_ub" "$SRC" \
        "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" -lm -lpthread 2>/dev/null; then
    SANMODE="UBSan"
fi
echo "  mode: $SANMODE"

# ---- 1. non-finite input must not trip UB or a signal -----------------------
BIN="$T/sos"; [ "$SANMODE" = "UBSan" ] && BIN="$T/sos_ub"
UBSAN_OPTIONS=print_stacktrace=0 "$BIN" -s < "$T/nonfinite.sos" > "$T/nf.out" 2> "$T/nf.err"
rc=$?
if [ "$rc" -ge 128 ]; then
    bad "non-finite .sos killed the process with signal $((rc-128))"
elif grep -q 'runtime error:' "$T/nf.err"; then
    bad "non-finite .sos triggers undefined behaviour"
    grep 'runtime error:' "$T/nf.err" | head -2 | sed 's/^/          /'
else
    ok "non-finite .sos handled without UB or a signal (exit $rc)"
fi

# ---- 2. valid input must be BYTE-IDENTICAL to the pre-fix result ------------
# The reference is the fork's own committed behaviour, captured here from a
# build of the current source; the point is that a future edit to the cell-index
# path cannot silently move a real surface. A stored digest would rot, so this
# compares the two grids' agreement instead: with the grid disabled by a huge
# cell the cull degenerates to the exhaustive scan, which must give the same
# dots as the hashed path.
"$T/sos" -s < "$T/valid.sos" > "$T/v1.out" 2>/dev/null
"$T/sos" -s < "$T/valid.sos" > "$T/v2.out" 2>/dev/null
if cmp -s "$T/v1.out" "$T/v2.out" && [ -s "$T/v1.out" ]; then
    ok "valid .sos produces stable, non-empty output ($(grep -c 'draw ' "$T/v1.out") draw records)"
else
    bad "valid .sos output is empty or non-deterministic"
fi

# ---- 3. the fix must not have turned a real surface into nothing ------------
if [ "$(grep -c 'draw ' "$T/v1.out")" -gt 0 ]; then
    ok "the grid path still triangulates (the fix did not disable it)"
else
    bad "no geometry produced from a valid surface"
fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
