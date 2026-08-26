#!/bin/sh
# HCAPEN cutoff-cache regression test.
#
# Builds test_hcapen_cache.f twice - against the patched hcapen_fast.f and
# against stock hcapen.f - and requires identical output. Both defects live in
# SAVE'd state across calls, so they are invisible to any single-shot check;
# verify.sh's whole-pipeline fixtures did not catch either one.
#
# Needs a hole2/src checkout for stock hcapen.f + its two dependencies. Skips
# cleanly (exit 0) when none is present, the same way the hydration parity test
# does for its gitignored inputs - but says so loudly rather than silently.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${1:-}"
if [ -z "$SRC" ]; then
  for c in "$HERE/stock_build/hole2/src" "$HERE/sphproc_build/hole2/src" "$HOME/hole2/src"; do
    [ -f "$c/hcapen.f" ] && SRC="$c" && break
  done
fi
if [ -z "$SRC" ] || [ ! -f "$SRC/hcapen.f" ]; then
  echo "hcapen_cache: SKIP - no hole2/src checkout with hcapen.f found."
  echo "              pass one as \$1 to run this test."
  exit 0
fi

FF="-fd-lines-as-comments -fbackslash -std=legacy -O2"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

gfortran $FF -o "$TMP/fast"  "$HERE/test_hcapen_cache.f" \
    "$HERE/hcapen_fast.f" "$SRC/holeen.f" "$SRC/ut_vector.f"
gfortran $FF -o "$TMP/stock" "$HERE/test_hcapen_cache.f" \
    "$SRC/hcapen.f" "$SRC/holeen.f" "$SRC/ut_vector.f"

"$TMP/fast"  > "$TMP/fast.txt"
"$TMP/stock" > "$TMP/stock.txt"

pass=0; fail=0
# Fortran's leading output space and column padding are not the thing under
# test - normalise whitespace on BOTH sides before comparing.
norm() { awk '{$1=$1};1'; }
norm < "$TMP/fast.txt"  > "$TMP/fast.n"
norm < "$TMP/stock.txt" > "$TMP/stock.n"
while read -r line; do
  tag="$(echo "$line" | awk '{print $1}')"
  want="$(grep "^$tag " "$TMP/stock.n" || true)"
  if [ "$line" = "$want" ]; then
    pass=$((pass+1)); echo "  PASS  $tag matches stock: $line"
  else
    fail=$((fail+1))
    echo "  FAIL  $tag"
    echo "        fast : $line"
    echo "        stock: $want"
  fi
done < "$TMP/fast.n"

# The cases are only meaningful if they actually exercise the defects: CASE1
# must find the LAST atom (past the 30000 cache cap) and CASE2B must find the
# moved interior atom. A driver that silently stopped doing that would pass
# the comparison above while testing nothing.
grep -q "CASE1 caprad= *1.000 iat1= *30001" "$TMP/stock.n" \
  || { echo "  FAIL  CASE1 no longer exercises the overflow path"; fail=$((fail+1)); }
grep -q "CASE2B caprad= *1.000 iat1= *2" "$TMP/stock.n" \
  || { echo "  FAIL  CASE2B no longer exercises the stale-cache path"; fail=$((fail+1)); }

echo "hcapen_cache: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
