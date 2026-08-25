#!/bin/sh
# Wrapper for hole_fast_coord.tcl - the packed coordinate record vs the PDB
# round trip. Never SKIPs at group level: with no patched binary it still
# asserts the gates that keep stock HOLE on the PDB path.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "hole-fast-coord: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" VMDHOLE_TEST_HOLE="${VMDHOLE_TEST_HOLE:-}" "$VMD" -dispdev text -e "$DIR/hole_fast_coord.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP|\.\.\.\.)'
res=$(echo "$out" | grep '^FASTCOORD-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
