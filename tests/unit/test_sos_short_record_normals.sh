#!/bin/sh
# A short .sos point record must read as "no normal", never as "the previous
# record's normal".
#
# THE DEFECT THIS GUARDS (verified to go red on the pre-fix source):
#   read_cord kept num5..num7 from the PREVIOUS sscanf when a record supplied
#   only 4 fields, on the strength of a comment claiming the per-vertex normals
#   are unconditionally recomputed later. They are not: vertex_normals()
#   returns before its assignment loop, so the file's normals flow through to
#   reorder_triangle(), which sums the three vertex normals and swaps two
#   corners when the dot product is negative. A short record therefore
#   inherited its predecessor's normal - and on the very first record, an
#   indeterminate stack value.
#
#   The fix zeroes num5..num7 before each sscanf, so a short record means
#   "no normal contribution". The invariant asserted here is exactly that: a
#   file with a short record must produce byte-identical output to the same
#   file with an explicit "0 0 0" normal on that record, even when the record
#   before it carries a normal that would flip triangle orientation if
#   inherited.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
SRC="$ROOT/native/sos_triangle_fast.c"
VOR="$ROOT/native/voronoi"
CC="${CC:-cc}"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "sos-short-record-normals: $SRC"
[ -f "$SRC" ] || { echo "SKIP: no sos_triangle_fast.c at $SRC"; exit 0; }
command -v "$CC" >/dev/null 2>&1 || { echo "SKIP: no C compiler ($CC)"; exit 0; }

T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

if ! "$CC" -O1 -o "$T/sos" "$SRC" "$VOR/vor_predicates.c" "$VOR/vor_delaunay.c" \
        -lm -lpthread 2>"$T/build.log"; then
    echo "SKIP: could not build sos_triangle_fast"; sed 's/^/    /' "$T/build.log" | head -5; exit 0
fi

# 60-dot synthetic surface, as in test_sos_nonfinite_input.sh. Record 30
# carries a strongly negative normal; record 31 is the probe. In "short" form
# it has 4 fields, in "zero" form an explicit "0 0 0" - the two files must be
# indistinguishable to the triangulator.
gen() { # $1 = short|zero for record 31, $2 = short|zero for record 0,
        # $3 = neg|zero for record 30's normal (default neg)
    awk -v m31="$1" -v r0="$2" -v n30="${3:-neg}" 'BEGIN{
        for (i = 0; i < 60; i++) {
            line = sprintf("4 %g %g %g", i*0.5, (i%7)*0.25, (i%5)*0.75)
            if      (i == 0  && r0  == "short") print line
            else if (i == 30 && n30 == "neg")   print line " -1 -1 -1"
            else if (i == 31 && m31 == "short") print line
            else                                print line " 0 0 0"
        }
    }'
}
gen short zero > "$T/mid_short.sos"
gen zero  zero > "$T/mid_zero.sos"
gen zero short > "$T/first_short.sos"
gen zero  zero > "$T/first_zero.sos"

run() { "$T/sos" -s < "$1" 2>/dev/null; }

run "$T/mid_short.sos" > "$T/a.out"
run "$T/mid_zero.sos"  > "$T/b.out"
ndraw=$(grep -c 'draw ' "$T/a.out")
if [ "$ndraw" -gt 0 ]; then
    ok "a surface containing a short record still triangulates ($ndraw draw records)"
else
    bad "short-record surface produced no geometry"
fi
if cmp -s "$T/a.out" "$T/b.out"; then
    ok "mid-file short record == explicit zero normal (no inheritance from record 30)"
else
    bad "mid-file short record differs from explicit zero: the previous normal leaked in"
fi

run "$T/first_short.sos" > "$T/c.out"
run "$T/first_zero.sos"  > "$T/d.out"
if cmp -s "$T/c.out" "$T/d.out"; then
    ok "first-record short record == explicit zero normal (no indeterminate value)"
else
    bad "first-record short record differs from explicit zero"
fi

# The equality checks above are one-sided: a build that zeroed EVERY normal
# (dropping file normals entirely) would pass them. Record 30 carrying
# -1 -1 -1 vs 0 0 0 must therefore change the output - normals really do
# reach reorder_triangle's vertex-order decision.
gen zero zero zero > "$T/all_zero.sos"
run "$T/all_zero.sos" > "$T/e.out"
if cmp -s "$T/b.out" "$T/e.out"; then
    bad "a nonzero file normal has no effect on output (normals dropped entirely?)"
else
    ok "a nonzero file normal changes the output (normals are consumed, not discarded)"
fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
