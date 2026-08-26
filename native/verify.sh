#!/bin/sh
# Reproduce the "byte-for-byte identical output" claim.
#
# Builds the unmodified upstream sos_triangle (../native/upstream) and the
# optimised sos_triangle_fast from source, runs both on every committed .sos
# test surface, and diffs the output. Needs only a C compiler -- no installed
# HOLE required.
#
# Part B additionally runs HOLE's own example structures through the full
# hole -> sph_process -> sos_triangle pipeline, IF the HOLE binaries are found
# (set HOLE_EXE to their directory; defaults to ~/hole2/exe).
#
# Usage:  ./verify.sh
set -e
cd "$(dirname "$0")"
CC="${CC:-cc}"; CFLAGS="${CFLAGS:--O2}"
# Stock reference: the pristine tree's binaries when they exist. ~/hole2/exe is
# where the ACCELERATED build gets installed, so defaulting to it made Part E
# compare the accelerated hole against itself once it was installed.
if [ -z "${HOLE_EXE:-}" ] && [ -x "$(dirname "$0")/stock_build/hole2/exe/hole" ]; then
  HOLE_EXE="$(cd "$(dirname "$0")" && pwd)/stock_build/hole2/exe"
fi
HOLE_EXE="${HOLE_EXE:-$HOME/hole2/exe}"
RAD="${RAD:-$HOME/hole2/rad/simple.rad}"
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT

# upstream, built two ways: as-shipped, and with only MAX_COORD raised so that
# surfaces which would overflow the 30000-polygon cap can still be compared.
# The unmodified upstream sos_triangle.c is not vendored (the repo is the
# plugin's code, not HOLE's): use a local copy if one exists, else fetch the
# file from the same pinned HOLE 2 revision build-vmdhole-optimized.sh clones.
UP="upstream/sos_triangle.c"
if [ ! -f "$UP" ]; then
    UP="$T/upstream_sos_triangle.c"
    curl -fsSL "https://raw.githubusercontent.com/osmart/hole2/master/src/sos_triangle.c" -o "$UP"         || { echo "SKIP: no local upstream/sos_triangle.c and the fetch failed - stock-vs-fast identity not checked"; exit 0; }
fi
$CC $CFLAGS -o "$T/stock" "$UP" -lm
sed 's/#define MAX_COORD 30000/#define MAX_COORD 200000/' "$UP" > "$T/big.c"
$CC $CFLAGS -o "$T/stockplus" "$T/big.c" -lm
$CC $CFLAGS -o "$T/fast" sos_triangle_fast.c -lm

pass=0; fail=0; overflow=0
diff_one () {   # $1 = .sos file ; reference picked automatically
  "$T/stock" -s < "$1" > "$T/o" 2>"$T/e" || true
  if grep -qi exceeded "$T/e"; then ref="$T/stockplus"; overflow=$((overflow+1)); else ref="$T/stock"; fi
  "$ref"     -s < "$1" > "$T/r" 2>/dev/null
  "$T/fast"  -s < "$1" > "$T/f" 2>/dev/null
  if diff -q "$T/r" "$T/f" >/dev/null; then pass=$((pass+1));
  else fail=$((fail+1)); echo "  MISMATCH: $1"; fi
}

echo "== Part A: committed fixtures + any local hole_output surfaces =="
# benchmarks/fixtures/*.sos is the committed corpus. hole_output/ is generated
# dev data (git-ignored), so it's diffed too when present locally but absent on
# a fresh clone.
# benchmarks/fixtures/ is GITIGNORED (generated locally), so on a fresh clone
# this glob matches nothing and expands to the literal pattern. Feeding that to
# diff_one produced a raw shell "cannot open" error and, under `set -e`, aborted
# the whole script - taking Parts B-E down with it even when they had every
# input they needed. Skip and continue, the way Parts B/C/E already do.
if [ -e benchmarks/fixtures/d15.sos ] || ls benchmarks/fixtures/*.sos >/dev/null 2>&1; then
  for s in benchmarks/fixtures/*.sos; do [ -f "$s" ] && diff_one "$s"; done
else
  echo "  skipped: no local corpus in benchmarks/fixtures/ (it is gitignored)."
  echo "           generate one, or run make_fixtures.sh, to enable Part A."
fi
for s in $(find ../vmdhole/hole_output -name '*.sos' 2>/dev/null); do diff_one "$s"; done
echo "  identical: $pass   mismatched: $fail   (of which needed raised cap: $overflow)"

echo "== Part B: HOLE example structures (full pipeline) =="
if [ ! -x "$HOLE_EXE/hole" ] || [ ! -x "$HOLE_EXE/sph_process" ] || [ ! -f "$RAD" ]; then
  echo "  skipped: set HOLE_EXE (=$HOLE_EXE) and RAD (=$RAD) to enable."
else
  EXDIR="${EXDIR:-/tmp/hole2_src/examples}"
  if [ ! -d "$EXDIR" ]; then
    echo "  cloning examples..."; git clone --depth 1 https://github.com/osmart/hole2.git "$T/hole2" >/dev/null 2>&1 && EXDIR="$T/hole2/examples"
  fi
  bp=0; bf=0
  for d in 01_gramicidin_1grm 02_choleratoxin_1chb 03_maltoporin_1af6; do
    [ -d "$EXDIR/$d" ] || continue
    w=$(mktemp -d); cp "$EXDIR/$d"/*.pdb "$w"/ 2>/dev/null || true
    inp=$(ls "$EXDIR/$d"/*.inp | head -1)
    sed "s|radius .*|radius $RAD|; s|sphpdb .*|sphpdb hole_out.sph|" "$inp" > "$w/h.inp"
    ( cd "$w" && "$HOLE_EXE/hole" < h.inp >/dev/null 2>&1; \
      "$HOLE_EXE/sph_process" -sos -dotden 15 hole_out.sph hole.sos >/dev/null 2>&1 ) || true
    [ -s "$w/hole.sos" ] || { echo "  $d: no .sos produced (skipped)"; rm -rf "$w"; continue; }
    "$T/stockplus" -s < "$w/hole.sos" > "$w/r" 2>/dev/null
    "$T/fast"      -s < "$w/hole.sos" > "$w/f" 2>/dev/null
    if diff -q "$w/r" "$w/f" >/dev/null; then echo "  $d: IDENTICAL"; bp=$((bp+1)); else echo "  $d: MISMATCH"; bf=$((bf+1)); fi
    rm -rf "$w"
  done
  echo "  examples identical: $bp   mismatched: $bf"
fi

echo "== Part C: hydrophobicity colouring (compiled --hydro vs Tcl reference) =="
# Confirms sos_triangle_fast's --hydro output matches the plugin's built-in Tcl
# colouring (ported in hydro_reference.py). Needs python3 + sph_process and the
# local hole_output/ frames (.sph + input_frame.pdb). That tree is git-ignored,
# so this part self-skips on a fresh clone.
SPHP="${SPHP:-./build/sph_process}"
[ -x "$SPHP" ] || SPHP="$HOLE_EXE/sph_process"
if ! command -v python3 >/dev/null 2>&1 || [ ! -x "$SPHP" ]; then
  echo "  skipped: need python3 and sph_process (set SPHP or HOLE_EXE)."
else
  hp=0; hf=0; found=0
  for sph in $(find ../vmdhole/hole_output -name hole_out.sph 2>/dev/null); do
    pdb="$(dirname "$sph")/input_frame.pdb"
    [ -f "$pdb" ] || continue
    found=$((found+1))
    "$SPHP" -sos -dotden 15 "$sph" "$T/h.sos" >/dev/null 2>&1; rm -f "$T/h.sos.old"
    "$T/fast" -s < "$T/h.sos" > "$T/base.plot" 2>/dev/null
    for sch in kd ww; do
      python3 hydro_reference.py sidecar  "$sph" "$pdb" $sch > "$T/atoms.dat"
      python3 hydro_reference.py colorize "$sph" "$T/atoms.dat" "$T/base.plot" $sch > "$T/ref.txt"
      "$T/fast" -s --hydro-atoms "$T/atoms.dat" --hydro-sph "$sph" --hydro-scheme $sch \
        < "$T/h.sos" > "$T/hyd.plot" 2>/dev/null
      awk '/^draw color/{c=$3} /draw trinorm/{print c}' "$T/hyd.plot" > "$T/fork.txt"
      d=$(paste "$T/fork.txt" "$T/ref.txt" | awk '$1!=$2{n++} END{print n+0}')
      if [ "$d" -eq 0 ]; then hp=$((hp+1)); else hf=$((hf+1)); echo "  MISMATCH $sph [$sch]: $d triangles"; fi
    done
  done
  if [ "$found" -eq 0 ]; then
    echo "  skipped: no frames with input_frame.pdb present (it is gitignored)."
  else
    echo "  colour-identical: $hp   mismatched: $hf"
    [ "$hf" -eq 0 ] || { echo "RESULT: hydrophobicity mismatch(es)."; exit 1; }
  fi
fi

echo "== Part D: dot surface (--points) =="
# --points must emit exactly the set of unique vertices of the trinorm mesh.
# Dependency-free: compares the two outputs of the fast binary on committed .sos.
dp=0; df=0
vtx_tri () { awk '/draw trinorm/{for(i=1;i<=NF;i++)gsub(/[{}]/,"",$i); n=0;
  for(i=1;i<=NF;i++) if($i ~ /^-?[0-9.]+$/) c[n++]=$i;
  for(j=0;j+2<n && j<9;j+=3) printf "%.3f %.3f %.3f\n",c[j],c[j+1],c[j+2]}' "$1" | sort -u; }
vtx_pts () { awk '/draw point/{for(i=1;i<=NF;i++)gsub(/[{}]/,"",$i); n=0;
  for(i=1;i<=NF;i++) if($i ~ /^-?[0-9.]+$/) c[n++]=$i;
  printf "%.3f %.3f %.3f\n",c[0],c[1],c[2]}' "$1" | sort -u; }
for s in benchmarks/fixtures/*.sos; do
  "$T/fast" -s          < "$s" > "$T/tri" 2>/dev/null
  "$T/fast" -s --points < "$s" > "$T/pts" 2>/dev/null
  vtx_tri "$T/tri" > "$T/tv"; vtx_pts "$T/pts" > "$T/pv"
  if diff -q "$T/tv" "$T/pv" >/dev/null; then dp=$((dp+1)); else df=$((df+1)); echo "  MISMATCH: $s"; fi
done
echo "  vertex-set identical: $dp   mismatched: $df"
[ "$df" -eq 0 ] || { echo "RESULT: dots mismatch(es)."; exit 1; }

echo "== Part E: parallel CONNOLLY / fast CAPSULE (accelerated vs stock hole) =="
# The accelerated `hole` (build-vmdhole-optimized.sh output) parallelises CONNOLLY
# and speeds up CAPSULE. CONN's per-plane progress lines are MUTED in the
# parallel path (shorto 2), so the comparison is the DATA that matters - the .sph point
# cloud (byte-identical) and the profile/conductance numbers - not the chatty
# progress log. Needs a stock hole (HOLE_EXE) and the accelerated one (ACCEL,
# default ./build/hole); self-skips otherwise.
ACCEL="${ACCEL:-./build/hole}"
if [ ! -x "$HOLE_EXE/hole" ] || [ ! -x "$ACCEL" ] || [ ! -f "$RAD" ]; then
  echo "  skipped: need stock \$HOLE_EXE/hole, accelerated \$ACCEL ($ACCEL), and \$RAD."
elif ldd "$HOLE_EXE/hole" 2>/dev/null | grep -qi gomp || nm "$HOLE_EXE/hole" 2>/dev/null | grep -q omp_fn; then
  # A stock reference that is itself an OpenMP build compares the accelerated
  # hole against itself and cannot fail. Refuse outright rather than pass.
  echo "  FAIL: \$HOLE_EXE/hole ($HOLE_EXE/hole) is an OpenMP (patched) build, not stock - point HOLE_EXE at a pristine build (native/stock_build/hole2/exe)."
  exit 1
else
  # Every run below happens inside "cd $w" (a scratch dir), so a RELATIVE
  # ACCEL/RAD (e.g. the "./build/hole" default) would resolve against $w, not
  # this script's directory, and silently fail to execute - "$ACCEL" not
  # found, swallowed by the trailing "|| true", producing a false
  # "MISMATCH (sph=no data=no)" instead of a mismatch OR a skip. Absolutize
  # both before the loop so the invocation is correct regardless of cwd.
  ACCEL=$(cd "$(dirname "$ACCEL")" && pwd)/$(basename "$ACCEL")
  RAD=$(cd "$(dirname "$RAD")" && pwd)/$(basename "$RAD")
  EXDIR="${EXDIR:-$T/hole2/examples}"
  [ -d "$EXDIR" ] || { git clone --depth 1 https://github.com/osmart/hole2.git "$T/hole2" >/dev/null 2>&1 && EXDIR="$T/hole2/examples"; }
  ep=0; ef=0
  for d in 01_gramicidin_1grm 03_maltoporin_1af6; do
    [ -d "$EXDIR/$d" ] || continue
    pdb=$(ls "$EXDIR/$d"/*.pdb 2>/dev/null | head -1); [ -f "$pdb" ] || continue
    src_inp=$(ls "$EXDIR/$d"/*.inp | head -1)
    for method in conn capsule; do
      w=$(mktemp -d); cp "$pdb" "$w/m.pdb"; cp "$RAD" "$w/r.rad"
      # bare filenames only (HOLE's 80-char FORTRAN line limit); run happens in $w.
      { echo "coord m.pdb"; echo "radius r.rad"; echo "sphpdb out.sph"; echo "raseed 42"; echo "shorto 1";
        grep -iE '^(cpoint|cvect|sample|endrad)' "$src_inp" 2>/dev/null || true; echo "$method"; echo "stop"; } > "$w/h.inp"
      ( cd "$w" && "$HOLE_EXE/hole" < h.inp > s.out 2>/dev/null; cp out.sph s.sph 2>/dev/null; rm -f out.sph out.sph.old ) || true
      ( cd "$w" && OMP_NUM_THREADS=8 OMP_STACKSIZE=256M "$ACCEL" < h.inp > a.out 2>/dev/null; cp out.sph a.sph 2>/dev/null ) || true
      # A run that produced NOTHING is not a mismatch - it is no evidence either
      # way, and reporting it as a mismatch is worse than useless: it fails the
      # headline byte-identity claim for a reason that has nothing to do with the
      # accelerated build. The usual cause is HOLE refusing the structure, e.g.
      # 1AF6 contains MG and the shipped simple.rad defines no MG radius, so HOLE
      # prints "***ERROR*** Cannot find vdW radius" and STILL EXITS 0. Part B
      # already treats absent output as a skip; Part E now does too.
      if [ ! -f "$w/s.sph" ] && [ ! -f "$w/a.sph" ]; then
        why=$(grep -A2 '\*\*\*ERROR\*\*\*' "$w/s.out" 2>/dev/null | tail -1 | tr -s ' ' | sed 's/^ //')
        echo "  $d [$method]: skipped - neither build produced a .sph${why:+ ($why)}"
        rm -rf "$w"; continue
      fi
      # .sph byte-identical, and profile ("(sampled)") + TAG lines identical
      sph_ok=no; dat_ok=no
      if [ -f "$w/s.sph" ] && [ -f "$w/a.sph" ] && diff -q "$w/s.sph" "$w/a.sph" >/dev/null 2>&1; then sph_ok=yes; fi
      grep -E '\(sampled\)|TAG ' "$w/s.out" > "$w/s.dat" 2>/dev/null || true
      grep -E '\(sampled\)|TAG ' "$w/a.out" > "$w/a.dat" 2>/dev/null || true
      if diff -q "$w/s.dat" "$w/a.dat" >/dev/null 2>&1; then dat_ok=yes; fi
      if [ "$sph_ok" = yes ] && [ "$dat_ok" = yes ]; then
        echo "  $d [$method]: IDENTICAL (.sph + profile + conductance)"; ep=$((ep+1))
      else
        echo "  $d [$method]: MISMATCH (sph=$sph_ok data=$dat_ok)"; ef=$((ef+1))
      fi
      rm -rf "$w"
    done
  done
  echo "  connolly/capsule identical: $ep   mismatched: $ef"
  [ "$ef" -eq 0 ] || { echo "RESULT: accelerated-hole mismatch(es)."; exit 1; }
fi

[ "$fail" -eq 0 ] && echo "RESULT: all committed surfaces byte-for-byte identical." || { echo "RESULT: $fail mismatch(es)."; exit 1; }
