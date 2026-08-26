#!/bin/sh
# remove_shallow must be bounded by the graph, not by the raw MinDepth value.
#
# THE DEFECT THIS GUARDS (measured on the pre-fix tree, 1ERI):
#   remove_shallow ran `for (d = min_depth; d >= 0; d--)` around a full nt scan,
#   twice per build, with an nt-int malloc/free inside the body. Every level
#   deeper than any live tetrahedron matches nothing - both loops select on
#   depth[i] == d - but still costs the scan. P.min_depth is an unclamped
#   atoi(argv[5]), and the GUI Min-depth field has no upper bound either.
#     min_depth      100000 ->   2122 ms
#     min_depth     5000000 -> 104900 ms   (linear; INT_MAX extrapolates to ~12 h)
#   min_depth + 1, used to seed the safe set, is also signed overflow at INT_MAX.
#
# The fix walks the depth levels that actually occur instead of the numeric
# range, which drops only empty levels and so cannot change the result.
set -u
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
ROOT=$(CDPATH= cd -- "$HERE/../.." && pwd)
ENG="$ROOT/native/mole_tunnel_engine"
AT="$ROOT/vmdhole/tests/fixtures/mole_atoms_1eri.txt"

pass=0; fail=0
ok()  { pass=$((pass+1)); echo "  PASS  $1"; }
bad() { fail=$((fail+1)); echo "  FAIL  $1"; }

echo "mole-depth-bound: $ENG"
[ -x "$ENG" ] || { echo "SKIP: no mole_tunnel_engine (run native/build.sh)"; exit 0; }
[ -f "$AT" ]  || { echo "SKIP: no 1ERI atom fixture"; exit 0; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT INT TERM

ms_run() {  # $1 = min_depth, $2 = output file -> echoes elapsed ms
    s=$(date +%s%N)
    timeout 120 "$ENG" "$AT" "$2" 3.0 1.25 "$1" 5 0 0 >/dev/null 2>&1
    rc=$?
    e=$(date +%s%N)
    echo "$(( (e-s)/1000000 )) $rc"
}

# --- the run must not scale with the parameter -------------------------------
set -- $(ms_run 8 "$T/base.out");        t_base=$1;  rc_base=$2
set -- $(ms_run 5000000 "$T/huge.out");  t_huge=$1;  rc_huge=$2

if [ "$rc_base" -ne 0 ]; then
    bad "the baseline run (min_depth=8) did not succeed (exit $rc_base)"
elif [ "$t_huge" -gt $(( t_base + 5000 )) ]; then
    bad "min_depth=5000000 took ${t_huge}ms vs ${t_base}ms at min_depth=8 - still scaling with the parameter"
else
    ok "min_depth=5000000 costs no more than the baseline (${t_huge}ms vs ${t_base}ms)"
fi

# --- INT_MAX must terminate, and not trip the min_depth+1 overflow -----------
set -- $(ms_run 2147483647 "$T/max.out"); t_max=$1; rc_max=$2
if [ "$rc_max" -eq 0 ] && [ "$t_max" -le $(( t_base + 5000 )) ]; then
    ok "min_depth=INT_MAX terminates cleanly (${t_max}ms, exit 0)"
else
    bad "min_depth=INT_MAX: ${t_max}ms exit $rc_max"
fi

# --- and the science must not have moved -------------------------------------
# Above the deepest level present, every larger MinDepth selects the same set,
# so results must be identical to each other. This is an invariant of the
# algorithm, so it needs no golden file to compare against.
"$ENG" "$AT" "$T/a.out" 3.0 1.25 500     5 0 0 >/dev/null 2>&1
"$ENG" "$AT" "$T/b.out" 3.0 1.25 5000000 5 0 0 >/dev/null 2>&1
if cmp -s "$T/a.out" "$T/b.out"; then
    ok "MinDepth values above the graph depth agree with each other"
else
    bad "MinDepth 500 and 5000000 gave different results"
fi

# A normal run must still find tunnels - guards against bounding it to nothing.
"$ENG" "$AT" "$T/n.out" 3.0 1.25 8 5 0 0 >/dev/null 2>&1
nt=$(grep -c '^T ' "$T/n.out" 2>/dev/null || echo 0)
[ "$nt" -gt 0 ] && ok "a default run still produces tunnels ($nt on 1ERI)" \
                || bad "a default run produced no tunnels"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
