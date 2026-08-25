#!/bin/sh
# Does hole::capsule agree with the real binary beyond the single fixture
# (1GRM, CVECT along z) every number in vmdhole/hole_tcl/README.md's capsule section
# was measured on?
#
# Three configurations, each seed-pinned (raseed 1, deterministic - HOLE's
# own MC noise only shows up across DIFFERENT seeds, see
# hole-monte-carlo-run-to-run-noise in this project's memory):
#   1. 1GRM, cvect 0 0 1   - the already-validated baseline, kept here as a
#      regression net now that a capsule test exists at all.
#   2. 1GRM, cvect 1 1 1   - a NON-axis-aligned CVECT. "1 1 1" is exactly
#      representable in float32 pre-normalisation (FREDA rounds every card
#      number through float32 - see README's FREDA trap note); HOLE
#      normalises AFTER that rounding, in double precision, same as
#      hole::capsule (capsule.tcl:307-309).
#   3. fixtures/1BL8_protein_only.pdb - a SECOND structure (KcsA, 4 chains,
#      2820 atoms, ~2x 1GRM) independent of 1GRM. No -cpoint given on either
#      side, so both engines derive it via CGUESS - exercising that path on
#      a new structure too. mcstep is turned down to 150 (from the real
#      default 1000, matched explicitly on BOTH sides via the `mcstep` card
#      / `-mcstep` flag) purely to keep this test's wall time reasonable
#      (~2820 atoms x 1000 steps/slice in the Tcl interpreter is ~2.3
#      minutes here - see vmdhole/hole_tcl/README.md's "Next steps" performance
#      note - vs ~20s at mcstep 150); it is still a real, fully deterministic
#      comparison at that reduced step count, not a shortcut that skips the
#      search. endrad is 6.0 (vs 1GRM's 8.0) - large enough to clear KcsA's
#      much narrower filter (bottleneck ~0.5-1 A above endrad 6 in <2 minutes)
#      without walking the whole ~60 A "atomic length" of this structure.
#
# NOTE: fixtures/1BL8_protein_only.pdb has its HETATM records (3 K+, 1 HOH)
# stripped - the ORIGINAL 1BL8.pdb makes the real binary abort ("Cannot find
# vdW radius for atom: K") because simple.rad has no K entry and tsatr.f has
# NO element fallback of its own (unlike this port's hole::element_radius,
# which is unconditional and would silently give K a radius instead of
# erroring - see that proc's own comment). See the fixture's REMARK header.
#
# This is a PRINT-PRECISION regression guard (the reference-build hole is not
# instrumented, so it can only be compared via its own %f-formatted stdout
# table) - it is NOT the evidence for the bit-exact claim in
# vmdhole/hole_tcl/README.md's "Root cause 3" section, which rests on a full E24.16
# trace of an instrumented Fortran tree (scratch-built, not checked in here -
# see that section for the reproducing recipe). What this test DOES catch:
# any future edit that regresses capsule.tcl's agreement with the real
# binary at the ~1e-4 A level, on a structure/CVECT this project's other
# capsule tests never exercise.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$DIR/../../.." && pwd)
REFBIN="$ROOT/vmdhole/hole_tcl/reference_bin"
RAD="$ROOT/vmdhole/tests/fixtures/simple.rad"
[ -f "$RAD" ] || RAD="$ROOT/native/stock_build/hole2/rad/simple.rad"

echo "=============================================================="
echo "capsule-generalization: hole::capsule vs the real binary beyond 1GRM/z"

[ -x "$REFBIN/hole" ] || {
  echo "  SKIP  no $REFBIN/hole - build it first:"
  echo "        cd native/stock_build/hole2/src && source ../source.apache && make ../exe/hole"
  echo "        then copy it into vmdhole/hole_tcl/reference_bin/ (local only, never committed - see vmdhole/hole_tcl/README.md)"
  exit 0
}
[ -f "$RAD" ] || { echo "  SKIP  no simple.rad found at $RAD"; exit 0; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
cp "$RAD" "$TMP/s.rad"

fail=0

# Compares the reference-build hole's own printed "(sampled)" rows (the capsule
# search's actual stored slices - NOT hcapgr.f's interpolated "(mid-point)"
# rows, which this port does not implement, see capsule.tcl's STATUS block)
# against hole.tcl -method capsule's -csv, index-for-index after sorting
# both by axial coordinate (slice counts must match exactly - they always
# have, once the fixes in capsule.tcl/hole.tcl landed - so this is a
# STRICTER pairing than profile_vs_reference.sh's nearest-t fuzzy match).
run_case () {
  label=$1; pdb=$2; cvect=$3; endrad=$4; mcstep=$5; cpoint=$6

  ref_out="$TMP/${label}_ref.txt"
  tcl_csv="$TMP/${label}_tcl.csv"
  bn=$(basename "$pdb")
  cp "$pdb" "$TMP/$bn"

  {
    printf 'coord %s\n' "$bn"
    printf 'radius s.rad\n'
    printf 'cvect %s\n' "$cvect"
    printf 'sample 0.5\n'
    printf 'endrad %s\n' "$endrad"
    printf 'raseed 1\n'
    printf 'mcstep %s\n' "$mcstep"
    [ -n "$cpoint" ] && printf 'cpoint %s\n' "$cpoint"
    printf 'capsul\n'
  } > "$TMP/${label}.inp"
  ( cd "$TMP" && "$REFBIN/hole" < "${label}.inp" > "${label}_ref.txt" 2>&1 )

  awk '/cenxyz.cvec/{f=1;next} f && /\(sampled\)/ && NF>=2 && ($1+0==$1) && ($2+0==$2) {printf "%.4f %.4f\n",$1,$2}' \
    "$ref_out" | sort -n | uniq > "$TMP/${label}_ref.dat"
  nref=$(wc -l < "$TMP/${label}_ref.dat")
  if [ "$nref" -eq 0 ]; then
    echo "  FAIL  $label: no reference (sampled) rows parsed"; tail -15 "$ref_out"; fail=1; return
  fi

  if [ -n "$cpoint" ]; then
    ( cd "$TMP" && tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb "$bn" -rad s.rad -method capsule \
        -cvect "$cvect" -sample 0.5 -endrad "$endrad" -seed 1 -mcstep "$mcstep" -cpoint "$cpoint" \
        -csv "${label}_tcl.csv" > "${label}_tcl.log" 2>&1 )
  else
    ( cd "$TMP" && tclsh "$ROOT/vmdhole/hole_tcl/hole.tcl" -pdb "$bn" -rad s.rad -method capsule \
        -cvect "$cvect" -sample 0.5 -endrad "$endrad" -seed 1 -mcstep "$mcstep" \
        -csv "${label}_tcl.csv" > "${label}_tcl.log" 2>&1 )
  fi
  [ -s "$tcl_csv" ] \
    || { echo "  FAIL  $label: hole.tcl errored"; cat "$TMP/${label}_tcl.log"; fail=1; return; }
  awk -F, 'NR>1{printf "%.4f %.4f\n",$1,$2}' "$tcl_csv" | sort -n > "$TMP/${label}_tcl.dat"
  ntcl=$(wc -l < "$TMP/${label}_tcl.dat")

  if [ "$nref" -ne "$ntcl" ]; then
    echo "  FAIL  $label: slice count mismatch, reference $nref vs hole.tcl $ntcl"
    fail=1; return
  fi

  paste "$TMP/${label}_ref.dat" "$TMP/${label}_tcl.dat" | awk -v label="$label" '
    { dz=$1-$3; if(dz<0)dz=-dz; dr=$2-$4; if(dr<0)dr=-dr
      if (dr>max) { max=dr; maxz=$1 }
      sum+=dr*dr; n++ }
    END {
      rms = sqrt(sum/n)
      printf "  %-22s %d slices, eff.rad RMS %.5f A, max %.5f A (at %.3f)\n", label, n, rms, max, maxz
      if (max > 0.0006) { print "  FAIL  " label ": max diff exceeds the 0.0006 A print-precision budget"; exit 1 }
    }' || fail=1
}

PDB1GRM="$ROOT/vmdhole/1GRM.pdb"
PDB1BL8="$DIR/fixtures/1BL8_protein_only.pdb"

if [ -f "$PDB1GRM" ]; then
  run_case "1grm_axis"    "$PDB1GRM" "0 0 1" 8.0 1000 "0 0 4"
  run_case "1grm_nonaxis" "$PDB1GRM" "1 1 1" 8.0 1000 "0 0 4"
else
  echo "  SKIP  no 1GRM.pdb found"
fi

if [ -f "$PDB1BL8" ]; then
  run_case "1bl8_2ndstruct" "$PDB1BL8" "0 0 1" 6.0 150 ""
else
  echo "  FAIL  fixture missing: $PDB1BL8"; fail=1
fi

[ "$fail" -eq 0 ] && echo "  PASS  all capsule-generalization cases agree with the real binary"
exit "$fail"
