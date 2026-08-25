#!/bin/sh
# Wrapper for capsule_incomplete.tcl - HOLE's own "RUN MAY BE INCOMPLETE"
# warning must reach the user. Needs the reference binary (the pure-Tcl
# fallback engine does not emit this warning); skips without it.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "capsule-incomplete: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" "$VMD" -dispdev text -e "$DIR/capsule_incomplete.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP|\.\.\.\.)'
res=$(echo "$out" | grep '^INCOMPLETE-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -5 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
