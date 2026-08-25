#!/bin/sh
# Acceleration-parity check: the accelerated binaries must BE accelerated, and
# must still produce byte-identical output.
#
# These are two different failures with two different fixes, so they are asserted
# separately:
#   * "accelerated"   - a build that silently lost its parallelism. This has now
#     happened three times (the R-15 class). The last one was invisible because
#     apply_patches.py's Makefile guard was all-or-nothing: a tree patched before
#     sphqpu_par was added could never receive it, so sph_process built SERIAL with
#     perfectly correct output. Hence the Makefile-invariant check below - it
#     catches that at PATCH time rather than after a build and a stopwatch.
#   * "identical"     - the whole justification for shipping the fork at all.
#
# Usage:  ./test_accel_parity.sh
#   EXE=<dir>    accelerated binaries      (default ~/hole2/exe)
#   STOCK=<file> stock sph_process to diff against (optional; skipped if absent)
#   SRC=<dir>    patched hole2/src tree to check Makefile invariants (optional)
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
EXE="${EXE:-$HOME/hole2/exe}"
# No default: a patched hole2/src tree is machine-local. Set SRC=<dir> to run
# the Makefile-invariant checks; unset, they skip cleanly.
SRC="${SRC:-}"
STOCK="${STOCK:-$DIR/../../native/stock_build/hole2/exe/sph_process}"
FIX="${FIX:-$DIR/fixtures}"
pass=0; fail=0; skip=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }
skp() { skip=$((skip+1)); echo "  SKIP  $1"; }

echo "accel-parity: EXE=$EXE"

# --- 1. Makefile invariants (catches the R-15 class at patch time) -------------
if [ -f "$SRC/Makefile" ]; then
  miss=""
  for o in coarea_fast.o hcapen_fast.o holcal_par.o holeen_par.o sphqpu_par.o; do
    grep -q "^$o" "$SRC/Makefile" || miss="$miss $o"
  done
  [ -z "$miss" ] && ok "all object swaps present in Makefile FILES" \
                 || bad "Makefile FILES missing object swap(s):$miss"
  miss=""
  for f in holcal_par holeen_par concal coarea_fast sphqpu_par; do
    grep -q "^$f\.o:" "$SRC/Makefile" || miss="$miss $f"
  done
  [ -z "$miss" ] && ok "all per-file -fopenmp rules present" \
                 || bad "Makefile missing per-file -fopenmp rule(s):$miss"
else
  skp "no patched src tree at $SRC (Makefile invariants)"
fi

# --- 2. the shipped binaries are actually OpenMP builds ------------------------
for b in hole sph_process; do
  if [ ! -x "$EXE/$b" ]; then skp "$b not found in $EXE"; continue; fi
  if ldd "$EXE/$b" 2>/dev/null | grep -qi gomp; then ok "$b is an OpenMP build"
  else bad "$b is NOT linked against libgomp - acceleration lost"; fi
done

# --- 3. sph_process actually scales (indicative, not a hard perf gate) ---------
BIG="$FIX/big.sph"
if [ -x "$EXE/sph_process" ] && [ -f "$BIG" ]; then
  T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
  ms() { s=$(date +%s%N); OMP_NUM_THREADS=$1 OMP_STACKSIZE=256M "$EXE/sph_process" \
           -sos -dotden 20 -color "$BIG" "$T/o$1.sos" >/dev/null 2>&1; e=$(date +%s%N)
         echo $(( (e-s)/1000000 )); }
  t1=$(ms 1); t8=$(ms 8)
  pct=$(( t8 * 100 / (t1>0 ? t1 : 1) ))
  if [ "$pct" -lt 70 ]; then ok "sph_process scales: t8 ${t8}ms = ${pct}% of t1 ${t1}ms"
  else bad "sph_process does NOT scale: t8 ${t8}ms = ${pct}% of t1 ${t1}ms (expect <70%)"; fi
else
  skp "scaling check (need $EXE/sph_process and $BIG)"
fi

# --- 4. byte-identical vs stock, across densities and thread counts ------------
# thread count 3 is deliberate: it does not divide the sphere count evenly, so a
# chunking bug in the parallel cull / ordered write shows up here and nowhere else.
if [ -x "$EXE/sph_process" ] && [ -x "$STOCK" ] && [ -f "$BIG" ]; then
  T2=$(mktemp -d); bad_n=0
  for d in 10 15 20 25; do
    OMP_NUM_THREADS=1 "$STOCK" -sos -dotden $d -color "$BIG" "$T2/ref.sos" >/dev/null 2>&1
    for t in 1 3 8; do
      OMP_NUM_THREADS=$t OMP_STACKSIZE=256M "$EXE/sph_process" \
        -sos -dotden $d -color "$BIG" "$T2/new.sos" >/dev/null 2>&1
      cmp -s "$T2/ref.sos" "$T2/new.sos" || { bad_n=$((bad_n+1)); echo "        differs: dotden=$d threads=$t"; }
    done
  done
  rm -rf "$T2"
  [ "$bad_n" -eq 0 ] && ok "sph_process byte-identical to stock (4 densities x threads 1,3,8)" \
                     || bad "$bad_n configuration(s) differ from stock"
else
  skp "byte-identity vs stock (need a stock sph_process at $STOCK)"
fi

# --- 5. CLI dispatch chain: a first-chain flag must not strand later flags ----
# Regression for the missing `else` that used to split main()'s flag-dispatch
# if/else-if chain in two (see the "BUG FIX" comment next to --tunnel-dist in
# sos_triangle_fast.c). Any flag matched in the FIRST half that does not itself
# return() (--asym-rays, --points, --hydro-atoms, ...) fell through into the
# SECOND half's own catch-all (help(); return(0)) on the very same argv[1], so
# every combined invocation silently printed the help banner and exited 0
# instead of doing the requested work - --asym-rays/--hydro-sph/--hydro3d-atoms
# combos were all affected, not just --asym-rays. test_accel_parity is the
# right home: it already runs against the DEPLOYED $EXE binary, which is what
# the user actually invokes - a source-only fix would not have caught this.
if [ -x "$EXE/sos_triangle" ]; then
  # Setup flags only - a terminal action flag (--asym-ellipse, --asymmetry) IS the
  # later flag this guards, so probing one with --hole-features proves nothing.
  for probe in "--recolor /dev/null" "--asym-rays 36" "--points" "--hydro-atoms /dev/null" \
               "--asym-threads 1"; do
    out=$("$EXE/sos_triangle" $probe --hole-features 2>/dev/null)
    case "$out" in
      hole_features:*) ok "sos_triangle $probe --hole-features reaches its handler" ;;
      *) bad "sos_triangle $probe --hole-features did NOT reach --hole-features (dispatch chain broken)" ;;
    esac
  done
else
  skp "CLI dispatch-chain check (need $EXE/sos_triangle)"
fi

echo "  -> $pass passed, $fail failed, $skip skipped"
[ "$fail" -eq 0 ]
