#!/bin/sh
# Wrapper for ellipse_parity.tcl - the ellipse C accelerator vs its Tcl reference.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "ellipse-parity: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" "$VMD" -dispdev text -e "$DIR/ellipse_parity.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^SKIP:|^  (PASS|FAIL|SKIP|\.\.\.\.)'
res=$(echo "$out" | grep '^ELLIPSE-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
