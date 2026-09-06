#!/bin/sh
# Nelder-Mead search + fast Connolly through the plugin's run path, compared
# with HOLE's Monte Carlo search and HOLE's own conn on 1GRM.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "nm-engine: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" VMDHOLE_TEST_HOLE="${VMDHOLE_TEST_HOLE:-}" "$VMD" -dispdev text -e "$DIR/nm_engine_e2e.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP|\.\.\.\.)|^SKIP:'
res=$(echo "$out" | grep '^NM-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -8 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
