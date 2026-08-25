#!/bin/sh
# Wrapper for hole_tcl_pore_methods.tcl - CONNOLLY/CAPSULE through the pure-Tcl
# fallback, compared against the reference binary. Skips without it.
# Slow by nature (~25 s per method: the fallback is the fallback).
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "hole-tcl-pore-methods: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" "$VMD" -dispdev text -e "$DIR/hole_tcl_pore_methods.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP|\.\.\.\.)'
res=$(echo "$out" | grep '^POREMETHOD-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
