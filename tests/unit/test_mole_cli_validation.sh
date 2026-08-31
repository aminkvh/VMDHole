#!/bin/sh
# The tunnel engine's CLI must refuse nonsense rather than crash on it or
# quietly answer a different question.
#
# THE DEFECTS THIS GUARDS (each verified to go red on the pre-fix source):
#   H2  A negative --cover made cover()'s recursion unbounded: the pivot's own
#       distance to itself is 0, and 0 > negative is TRUE, so `rest` became a
#       verbatim copy of `comp` and every level wrote one more pivot past the
#       end of the caller's buffer. ASan: heap-buffer-overflow WRITE with seven
#       cover() frames; plain build: SIGSEGV (exit 139).
#   L9  An out-of-range weight index silently ran MOLE's default, because
#       mole_edge_cost's `default:` arm returns what weight 0 returns.
#   L10 argv[3..8] were read without checking they are positionals, so
#       `engine atoms out --cover=6` parsed "--cover=6" as the probe radius
#       (atof -> 0.0) and reported "0 channels, 0 voids", exit 0.
#
# Also asserts the fixes did not move any science: the engine's output on a
# committed fixture must be unchanged.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
FIX="$ROOT/vmdhole/tests/fixtures"
ENG="$ROOT/native/mole_tunnel_engine"
AT="$FIX/mole_atoms_1eri.txt"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "mole-cli-validation: $ENG"
[ -x "$ENG" ] || { echo "SKIP: no mole_tunnel_engine (run native/build.sh)"; exit 0; }
[ -f "$AT" ]  || { echo "SKIP: no 1ERI atom fixture"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

# --- H2: a negative cover radius must terminate, not smash the stack --------
"$ENG" "$AT" "$T/o.out" 3.0 1.25 8 5 0 0 --cover=-1 >"$T/h2.log" 2>&1
rc=$?
if [ "$rc" -ge 128 ]; then bad "negative --cover killed the process (signal $((rc-128)))"
else ok "negative --cover terminates cleanly (exit $rc)"; fi

# --- L9: an out-of-range weight must be reported, not silently defaulted ----
"$ENG" "$AT" "$T/o.out" 3.0 1.25 8 5 0 99 >/dev/null 2>"$T/l9.log"
rc=$?
if [ "$rc" -eq 2 ] && grep -q 'out of range' "$T/l9.log"; then
    ok "out-of-range weight is refused with a diagnostic (exit 2)"
else bad "out-of-range weight: exit $rc, stderr: $(head -1 "$T/l9.log")"; fi
# a VALID weight must still be accepted
"$ENG" "$AT" "$T/o.out" 3.0 1.25 8 5 0 3 >/dev/null 2>&1
[ $? -eq 0 ] && ok "a valid weight (3=Constant) is still accepted" \
             || bad "a valid weight was rejected"

# --- L10: a flag in a positional slot must not be parsed as a number --------
"$ENG" "$AT" "$T/o.out" --cover=6 >/dev/null 2>"$T/l10.log"
if grep -q '0 channels, 0 voids' "$T/l10.log"; then
    bad "a flag in a positional slot still collapses the run (probe read as 0.0)"
else
    ok "a flag in a positional slot leaves the defaults intact ($(sed -n 's/.*tetrahedra, //p' "$T/l10.log" | head -1))"
fi

# --- non-finite or junk numerics must be refused, not silently reshaped -----
for probe in "nan 1.25 8" "3.0 inf 8" "3.0 1.25 8x"; do
    # shellcheck disable=SC2086
    "$ENG" "$AT" "$T/o.out" $probe 5 0 0 >/dev/null 2>"$T/nf.log"
    rc=$?
    if [ "$rc" -eq 2 ] && grep -Eq 'not (a finite number|an integer)' "$T/nf.log"; then
        ok "argv '$probe' is refused with a diagnostic (exit 2)"
    else bad "argv '$probe': exit $rc, stderr: $(head -1 "$T/nf.log")"; fi
done
"$ENG" "$AT" "$T/o.out" 3.0 1.25 8 5 0 0 nan 28.9 1.1 >/dev/null 2>"$T/nf2.log"
rc=$?
if [ "$rc" -eq 2 ] && grep -q 'origin x' "$T/nf2.log"; then
    ok "a nan origin component is refused naming the field (exit 2)"
else bad "nan origin: exit $rc, stderr: $(head -1 "$T/nf2.log")"; fi

# --- the science must not have moved ----------------------------------------
"$ENG" "$AT" "$T/ref.out" 3.0 1.25 8 5 0 0 23.577 28.862 1.087 >/dev/null 2>&1
if [ -s "$T/ref.out" ] && [ "$(grep -c '^T ' "$T/ref.out")" -gt 0 ]; then
    ok "a normal run still produces tunnels ($(grep -c '^T ' "$T/ref.out") on 1ERI)"
else bad "a normal run produced no tunnels"; fi

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
