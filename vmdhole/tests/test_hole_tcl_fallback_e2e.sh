#!/bin/sh
# Wrapper for hole_tcl_fallback_e2e.tcl - the fallback through run_analysis
# itself, under `vmd -dispdev text` (run_analysis has a documented headless
# batch path, so this needs no display).
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "hole-tcl-fallback-e2e: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" "$VMD" -dispdev text -e "$DIR/hole_tcl_fallback_e2e.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^SKIP:|^  (PASS|FAIL|SKIP|\.\.\.\.)'
res=$(echo "$out" | grep '^E2E-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
