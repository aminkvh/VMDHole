#!/bin/sh
# A coordinate that is not a finite number must stop the run, not become 0.
#
# THE DEFECT THIS GUARDS:
#   Two readers, two ways to lose a coordinate silently.
#     * mole_parse_double faithfully reproduces MOLE's NumberParser.
#       ParseDoubleFast, which has no NaN/Inf handling - it stops at the first
#       character it does not recognise, so "nan" parses as 0.0 and the atom
#       lands at the origin. The parser is deliberately NOT corrected: bit
#       fidelity with MOLE is what it exists for.
#     * the MOLE_ATOMS_EXACT path uses strtod, which DOES accept "nan"/"inf",
#       and the non-finite value then propagates into the Delaunay predicates.
#   Either way the engine exits 0 and reports a plausible tunnel set built from
#   a coordinate nobody supplied - the same outcome mole_read_atoms already
#   refuses for a truncated structure.
#
# The fix validates the FIELD and leaves the parser alone, so every well-formed
# file is bit-for-bit unchanged. This test checks both halves.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
ENG="$ROOT/native/mole_tunnel_engine"
AT="$ROOT/vmdhole/tests/fixtures/mole_atoms_1eri.txt"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "mole-coord-validation: $ENG"
[ -x "$ENG" ] || { echo "SKIP: no mole_tunnel_engine (run native/build.sh)"; exit 0; }
[ -f "$AT" ]  || { echo "SKIP: no 1ERI atom fixture"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# --- every well-formed fixture must still run, and be unchanged -------------
nfix=0; nrun=0
for f in "$ROOT"/vmdhole/tests/fixtures/mole_atoms_*.txt; do
    [ -f "$f" ] || continue
    nfix=$((nfix+1))
    if "$ENG" "$f" "$T/ok.out" 3.0 1.25 8 5 0 0 >/dev/null 2>&1; then nrun=$((nrun+1)); fi
done
if [ "$nfix" -gt 0 ] && [ "$nrun" -eq "$nfix" ]; then
    ok "all $nfix well-formed fixtures still run (validation rejects nothing valid)"
else
    bad "$((nfix-nrun)) of $nfix well-formed fixtures were rejected"
fi

# --- a non-finite or non-numeric coordinate must be refused -----------------
# Column 1 of the fixture is the x coordinate; replace it on one data row.
for badv in nan NaN inf -inf abc; do
    awk -v b="$badv" 'NR==3 { sub(/^[^ \t]+/, b) } { print }' "$AT" > "$T/bad.txt"
    if [ "$(head -3 "$T/bad.txt" | tail -1 | cut -d' ' -f1)" != "$badv" ]; then
        bad "could not build the $badv fixture - the substitution stopped matching"
        continue
    fi
    out=$("$ENG" "$T/bad.txt" "$T/bad.out" 3.0 1.25 8 5 0 0 2>&1); rc=$?
    if [ "$rc" -ne 0 ] && printf '%s' "$out" | grep -q 'not a finite number'; then
        ok "coordinate \"$badv\" is refused with a diagnostic (exit $rc)"
    else
        bad "coordinate \"$badv\": exit $rc, no diagnostic (first line: $(printf '%s' "$out" | head -1 | cut -c1-60))"
    fi
done

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
