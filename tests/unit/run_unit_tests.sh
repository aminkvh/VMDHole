#!/bin/sh
# Run every unit test in this directory and report a total.
#
# These are deliberately SMALL and SELF-CONTAINED: no VMD, no trajectory data,
# no gitignored fixture corpus, no network. Each builds or generates whatever it
# needs into a temp dir and cleans up. That is what makes them safe to run on
# every commit, unlike the main suite whose heavier groups need a licensed VMD
# and locally built HOLE binaries.
#
# Exit 0 = all passed. A test that genuinely cannot run prints "SKIP:" at
# column 0 and exits 0 - the same contract vmdhole/tests/run_tests.sh uses.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SELF=$(basename -- "$0")
pass=0; fail=0; skip=0
for t in "$DIR"/test_*.sh; do
    [ -e "$t" ] || continue
    name=$(basename -- "$t")
    [ "$name" = "$SELF" ] && continue
    echo "=============================================================="
    out=$("$t" 2>&1); rc=$?
    echo "$out"
    if [ $rc -ne 0 ]; then fail=$((fail+1)); echo "  >>> $name: FAIL (exit $rc)"
    elif printf '%s\n' "$out" | grep -q '^SKIP:'; then skip=$((skip+1))
    else pass=$((pass+1)); fi
done
echo "=============================================================="
echo "unit: $pass passed, $fail failed, $skip skipped"
[ "$fail" -eq 0 ]
