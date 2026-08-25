#!/bin/sh
# Headless smoke test: the plugin must load, parse, and import under `vmd -dispdev text`,
# which has NO Tk. Wraps headless_smoke.tcl and turns its counters into an exit status.
# Usage:  ./test_headless_smoke.sh      (VMD=<path> to override the vmd binary)
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "headless-smoke: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" "$VMD" -dispdev text -e "$DIR/headless_smoke.tcl" 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP)|^        ' 
res=$(echo "$out" | grep '^SMOKE-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end (plugin load or VMD failure)"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
