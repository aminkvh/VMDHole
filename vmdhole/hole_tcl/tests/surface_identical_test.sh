#!/bin/sh
# Does the Tcl 3D-surface pipeline (sph_process.tcl + sos_triangle.tcl)
# produce BYTE-IDENTICAL .sos and .plot files to the real binaries?
#
# Uses vmdhole/hole_tcl/reference_bin/{hole,sph_process,sos_triangle} (local, never committed) - a build from
# native/stock_build/hole2/src's pristine, unpatched Makefile
# (`source ../source.apache && make`), NOT ~/hole2/exe. ~/hole2/exe/hole was
# found to be linked against libgomp with a holcal_._omp_fn.0 symbol: it is a
# build of native/connolly_patches's experimental parallel-CONNOLLY
# patch set, not stock HOLE, and SEGFAULTS on CONNOLLY+RASEED (reproduced by
# applying those same patches to a clean tree and rebuilding). It must never
# be used as a fidelity reference again. See vmdhole/hole_tcl/README.md (Reference binaries section).
#
# Only sph_process/sos_triangle are exercised here (not `hole` itself): this
# isolates the surface pipeline from the annealing search's own Monte Carlo
# noise, matching how it was originally verified.
#
# Deliberately does NOT source hole.tcl or capsule.tcl - both are being
# actively edited elsewhere and sourcing a mid-edit file would make this test
# flaky by construction. sph_process.tcl/sos_triangle.tcl are sourced
# directly; sos_triangle.tcl only needs hole::_f32 from sph_process.tcl.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$DIR/../../.." && pwd)
REFBIN="$ROOT/vmdhole/hole_tcl/reference_bin"

echo "=============================================================="
echo "surface-identical: sph_process.tcl + sos_triangle.tcl vs the fresh binaries"

for b in hole sph_process sos_triangle; do
  [ -x "$REFBIN/$b" ] || {
    echo "  SKIP  no $REFBIN/$b - build it first:"
    echo "        cd native/stock_build/hole2/src && source ../source.apache && make ../exe/hole ../exe/sph_process ../exe/sos_triangle"
    echo "        then copy the three binaries into vmdhole/hole_tcl/reference_bin/ (local only - see vmdhole/hole_tcl/README.md)"
    exit 0
  }
done

PDB="$ROOT/vmdhole/1GRM.pdb"
[ -f "$PDB" ] || { echo "  SKIP  no 1GRM.pdb found"; exit 0; }
RAD="$ROOT/vmdhole/tests/fixtures/simple.rad"
[ -f "$RAD" ] || RAD="$ROOT/native/stock_build/hole2/rad/simple.rad"
[ -f "$RAD" ] || { echo "  SKIP  no simple.rad found"; exit 0; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
# Short filenames on purpose: HOLE's card fields are fixed width and silently
# truncate a long path, which reads as "cannot open radius file".
cp "$PDB" "$TMP/in.pdb"
cp "$RAD" "$TMP/s.rad"
cd "$TMP" || exit 1

# Generate the .sph with the real hole binary (raseed pinned - irrelevant to
# what follows, but keeps the run reproducible if inspected by hand).
cat > ctrl.inp <<EOF
coord in.pdb
radius s.rad
sphpdb out.sph
cvect 0.0 0.0 1.0
sample 0.5
endrad 8.0
raseed 1
EOF
"$REFBIN/hole" < ctrl.inp > hole_run.log 2>&1
[ -s out.sph ] || { echo "  FAIL  reference hole did not produce out.sph"; tail -20 hole_run.log; exit 1; }

# --- .sos: sph_process (real binary) vs hole::sph_process (Tcl) ---
"$REFBIN/sph_process" -sos -dotden 10 out.sph bin.sos > sph_process.log 2>&1
[ -s bin.sos ] || { echo "  FAIL  reference sph_process produced no output"; cat sph_process.log; exit 1; }

cat > run_sph.tcl <<EOF
source "$ROOT/vmdhole/hole_tcl/sph_process.tcl"
hole::sph_process out.sph tcl.sos 10
EOF
tclsh run_sph.tcl > tcl_sph.log 2>&1 || { echo "  FAIL  hole::sph_process errored"; cat tcl_sph.log; exit 1; }

if diff -q bin.sos tcl.sos > /dev/null 2>&1; then
  echo "  PASS  .sos byte-identical ($(wc -l < bin.sos) lines)"
else
  echo "  FAIL  .sos differs"
  diff bin.sos tcl.sos | head -10
  exit 1
fi

# --- .plot: sos_triangle -s (real binary) vs hole::sos_triangle (Tcl) ---
"$REFBIN/sos_triangle" -s < bin.sos > bin.plot 2> sos_triangle.log
[ -s bin.plot ] || { echo "  FAIL  reference sos_triangle produced no output"; cat sos_triangle.log; exit 1; }

cat > run_sos.tcl <<EOF
source "$ROOT/vmdhole/hole_tcl/sph_process.tcl"
source "$ROOT/vmdhole/hole_tcl/sos_triangle.tcl"
hole::sos_triangle bin.sos tcl.plot
EOF
tclsh run_sos.tcl > tcl_sos.log 2>&1 || { echo "  FAIL  hole::sos_triangle errored"; cat tcl_sos.log; exit 1; }

if diff -q bin.plot tcl.plot > /dev/null 2>&1; then
  echo "  PASS  .plot byte-identical ($(wc -l < bin.plot) lines, $(grep -c '^draw trinorm' bin.plot) triangles)"
else
  echo "  FAIL  .plot differs"
  diff bin.plot tcl.plot | head -10
  exit 1
fi

echo "  PASS  end-to-end .sph -> .sos -> .plot byte-identical against the fresh binaries"
exit 0
