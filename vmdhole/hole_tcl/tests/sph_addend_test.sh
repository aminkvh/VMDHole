#!/bin/sh
# Does hole::write_sph (and its capsule/connolly siblings) produce a
# BYTE-IDENTICAL .sph to the reference-build hole's own `sphpdb` output - discovery
# order, IREC numbering, the ADDEND terminal-sphere grid and its
# LAST-REC-END markers, and the duplicate record-0 rewrite at the +ve/-ve
# transition - not just a .sph that merely parses?
#
# Before this test existed, hole::write_sph wrote per-slice ATOM records
# SORTED by axial coordinate, numbered 0..N-1 positionally (not the real
# STRNOP/-STRNON IREC), and never wrote ADDEND's grid or LAST-REC-END at all
# - 82 lines vs the real binary's 421 on the fixture this project's other
# fidelity numbers use (1GRM, cpoint 0 0 0, cvect 0 0 1, sample 0.25, endrad
# 22.0, raseed 1). See hole::addend's and hole::write_sph's own headers in
# hole.tcl for the full derivation (addend.f, wpdbsp.f, holcal.f ~758-806).
#
# mcstep is turned down to 50 (from the real default 1000) purely to keep
# this test's wall time proportionate to the rest of run_all.sh - matched
# explicitly on BOTH sides via the `mcstep` card / `-mcstep` flag, so it is
# still a real, fully deterministic byte-comparison at that reduced step
# count (HOLE's RNG stream is bit-exact regardless of mcstep - see
# vmdhole/hole_tcl/README.md's RNG section) - endrad/sample/cpoint are chosen so both
# directions still store >=2 slices each at this step count (ADDEND's own
# LASCEN reset needs that - see write_sph's header). Connolly's own grid is
# widened to 0.6 (from the project's usual 0.2) for the same reason - its
# cost is governed by CONNR(2), not mcstep, and this test only needs the
# record FORMAT/ORDER to be exercised, not a fine-grained area estimate.
#
# Measured cost of this file alone (see the delegate report that added it):
# spherical ~1.5s, capsule ~2s, connolly ~8s, radius-abort check <1s wall -
# about 12s added to run_all.sh's existing ~75s.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$DIR/../../.." && pwd)
REFBIN="$ROOT/vmdhole/hole_tcl/reference_bin"

echo "=============================================================="
echo "sph-addend: hole::write_sph/write_capsule_sph/write_connolly_sph vs the real binary's sphpdb output"

[ -x "$REFBIN/hole" ] || {
  echo "  SKIP  no $REFBIN/hole - build it first:"
  echo "        cd native/stock_build/hole2/src && source ../source.apache && make ../exe/hole"
  echo "        then copy it into vmdhole/hole_tcl/reference_bin/ (local only - see vmdhole/hole_tcl/README.md)"
  exit 0
}

PDB="$ROOT/vmdhole/1GRM.pdb"
[ -f "$PDB" ] || { echo "  SKIP  no 1GRM.pdb found"; exit 0; }
RAD="$ROOT/vmdhole/tests/fixtures/simple.rad"
[ -f "$RAD" ] || RAD="$ROOT/native/stock_build/hole2/rad/simple.rad"
[ -f "$RAD" ] || { echo "  SKIP  no simple.rad found"; exit 0; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cp "$PDB" "$TMP/in.pdb"
cp "$RAD" "$TMP/s.rad"
cd "$TMP" || exit 1

fail=0

# Common cards: cpoint 0 0 4, cvect 0 0 1, sample 0.5, endrad 8.0, raseed 1,
# mcstep 50 - shared by all three methods below so one .inp template covers
# them (method-specific cards appended per case).
common () {
  printf 'coord in.pdb\nradius s.rad\ncpoint 0 0 4\ncvect 0 0 1\n'
  printf 'sample 0.5\nendrad 8.0\nraseed 1\nmcstep 50\n'
}

# --- spherical ---
{ common; printf 'sphpdb ref_sph.sph\n'; } > sph.inp
"$REFBIN/hole" < sph.inp > sph_ref.log 2>&1
[ -s ref_sph.sph ] || { echo "  FAIL  spherical: reference hole did not produce ref_sph.sph"; tail -15 sph_ref.log; fail=1; }
tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb in.pdb -rad s.rad -cpoint "0 0 4" -cvect "0 0 1" \
    -sample 0.5 -endrad 8.0 -seed 1 -mcstep 50 -sph tcl_sph.sph > sph_tcl.log 2>&1
if [ -s tcl_sph.sph ] && diff -q ref_sph.sph tcl_sph.sph > /dev/null 2>&1; then
  echo "  PASS  spherical .sph byte-identical ($(wc -l < ref_sph.sph) lines: $(grep -c '^ATOM' ref_sph.sph) ATOM, $(grep -c 'LAST-REC-END' ref_sph.sph) LAST-REC-END)"
else
  echo "  FAIL  spherical .sph differs"
  cat sph_tcl.log
  diff ref_sph.sph tcl_sph.sph 2>&1 | head -10
  fail=1
fi

# --- capsule ---
{ common; printf 'sphpdb ref_cap.sph\ncapsul\n'; } > cap.inp
"$REFBIN/hole" < cap.inp > cap_ref.log 2>&1
[ -s ref_cap.sph ] || { echo "  FAIL  capsule: reference hole did not produce ref_cap.sph"; tail -15 cap_ref.log; fail=1; }
tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb in.pdb -rad s.rad -method capsule -cpoint "0 0 4" -cvect "0 0 1" \
    -sample 0.5 -endrad 8.0 -seed 1 -mcstep 50 -sph tcl_cap.sph > cap_tcl.log 2>&1
if [ -s tcl_cap.sph ] && diff -q ref_cap.sph tcl_cap.sph > /dev/null 2>&1; then
  echo "  PASS  capsule .sph byte-identical ($(wc -l < ref_cap.sph) lines: $(grep -c QC1 ref_cap.sph) QC1/$(grep -c QC2 ref_cap.sph) QC2)"
else
  echo "  FAIL  capsule .sph differs"
  cat cap_tcl.log
  diff ref_cap.sph tcl_cap.sph 2>&1 | head -10
  fail=1
fi

# --- connolly --- (grid widened to 0.6, see header)
{ common; printf 'sphpdb ref_conn.sph\nconn 1.15 0.6\n'; } > conn.inp
"$REFBIN/hole" < conn.inp > conn_ref.log 2>&1
[ -s ref_conn.sph ] || { echo "  FAIL  connolly: reference hole did not produce ref_conn.sph"; tail -15 conn_ref.log; fail=1; }
tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb in.pdb -rad s.rad -method connolly -cpoint "0 0 4" -cvect "0 0 1" \
    -sample 0.5 -endrad 8.0 -seed 1 -mcstep 50 -probe 1.15 -grid 0.6 -sph tcl_conn.sph > conn_tcl.log 2>&1
if [ -s tcl_conn.sph ] && diff -q ref_conn.sph tcl_conn.sph > /dev/null 2>&1; then
  echo "  PASS  connolly .sph byte-identical ($(wc -l < ref_conn.sph) lines)"
else
  echo "  FAIL  connolly .sph differs"
  cat conn_tcl.log
  diff ref_conn.sph tcl_conn.sph 2>&1 | head -10
  fail=1
fi

# --- element_radius fallback: the real binary ABORTS on an atom no VDWR
# rule matches (tsatr.f, no fallback of its own); hole::radius_for must now
# do the same instead of silently substituting a carbon radius - see
# hole::radius_for/hole::element_radius's own comments in hole.tcl.
K_PDB="$DIR/fixtures/1BL8_protein_only.pdb"
if [ -f "$K_PDB" ]; then
  cat > kpdb.pdb <<EOF
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  K   K   A   2       5.000   0.000   0.000  1.00  0.00           K
EOF
  if tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb kpdb.pdb -rad s.rad -cpoint "0 0 0" \
      -cvect "0 0 1" -sample 1.0 -endrad 5.0 -mcstep 5 -csv kpdb.csv \
      > kpdb.log 2>&1; then
    echo "  FAIL  radius_for: K+ with simple.rad should have aborted like the real binary (tsatr.f has no element fallback) but exited 0"
    fail=1
  elif grep -q "Cannot find vdW radius for atom" kpdb.log; then
    echo "  PASS  radius_for aborts on an unmatched VDWR atom, matching tsatr.f's own message"
  else
    echo "  FAIL  radius_for: errored but not with tsatr.f's own message"
    cat kpdb.log
    fail=1
  fi
else
  echo "  SKIP  no $K_PDB - radius_for abort check needs it"
fi

[ "$fail" -eq 0 ] && echo "  PASS  all .sph writers agree with the real binary's sphpdb output"
exit "$fail"
