#!/bin/sh
# The marching-cubes spherical mesher (mesh_csg) behind the surface display:
# gating, naming, persistent server, accuracy, and the sos_triangle fallback.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
VMD="${VMD:-vmd}"
command -v "$VMD" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "mesh-csg-engine: $VMD -dispdev text"
out=$(VMDHOLE_TEST_DIR="$DIR" VMDHOLE_TEST_HOLE="${VMDHOLE_TEST_HOLE:-}" "$VMD" -dispdev text -e "$DIR/mesh_csg_engine_e2e.tcl" < /dev/null 2>&1)
echo "$out" | grep -E '^  (PASS|FAIL|SKIP|\.\.\.\.)|^SKIP:'
res=$(echo "$out" | grep '^MC-RESULT')
[ -n "$res" ] || { echo "  FAIL  test did not reach the end"; \
                   echo "$out" | tail -8 | sed 's/^/        /'; exit 1; }
f=$(echo "$res" | sed 's/.*fail=//')
p=$(echo "$res" | sed 's/.*pass=\([0-9]*\).*/\1/')
echo "  -> $p passed, $f failed"
[ "$f" -eq 0 ]
