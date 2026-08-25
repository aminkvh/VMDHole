#!/bin/sh
# Run every test in this directory and report a total.
#
# Selects by *_test.sh, not by executable bit or a hardcoded list: this
# deliberately excludes profile_vs_reference.sh, which takes a required
# positional PDB argument (usage: profile_vs_reference.sh PDB ...) and is a
# harness invoked by hand, not a self-contained pass/fail test - running it
# with no args would just print its usage error and be counted as a bogus
# FAIL.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SELF=$(basename -- "$0")

pass=0
fail=0
skip=0

for t in "$DIR"/*_test.sh; do
  [ -e "$t" ] || continue
  name=$(basename -- "$t")
  [ "$name" = "$SELF" ] && continue
  out=$("$t" < /dev/null 2>&1)
  rc=$?
  echo "$out"
  if [ $rc -ne 0 ]; then
    fail=$((fail+1))
    echo "  >>> $name: FAIL (exit $rc)"
  elif printf '%s\n' "$out" | grep -q '^  SKIP'; then
    skip=$((skip+1))
  else
    pass=$((pass+1))
  fi
  echo
done

echo "=============================================================="
echo "run_all: $pass passed, $fail failed, $skip skipped"
[ $fail -eq 0 ]
