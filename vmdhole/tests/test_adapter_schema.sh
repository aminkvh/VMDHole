#!/bin/sh
# Runs adapter_schema.tcl under VMD. Needs Tcl 8.6, so it prefers vmd2 and
# SKIPS where the GUI-capable build is unavailable - the adapter itself is
# Tk-free, but the plugin declines to load without Tk 8.6.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
LOG="${TMPDIR:-/tmp}/vmdhole_adapter.$$.log"
BIN="${VMD_BIN:-}"
if [ -z "$BIN" ]; then
    for c in vmd2 vmd; do command -v "$c" >/dev/null 2>&1 && { BIN=$c; break; }; done
fi
[ -n "$BIN" ] || { echo "SKIP: no vmd on PATH"; exit 0; }
echo "adapter-schema: $BIN"
ADAPTER_LOG="$LOG" ADAPTER_PDB="$DIR/fixtures/altloc_icode.pdb" ADAPTER_HET_PDB="$DIR/fixtures/hetatm_all.pdb" \
VMDHOLE_TCL="$DIR/../vmdhole.tcl" \
    timeout 180 "$BIN" -dispdev text -e "$DIR/adapter_schema.tcl" >/dev/null 2>&1
[ -s "$LOG" ] || { echo "SKIP: adapter check did not run under $BIN"; rm -f "$LOG"; exit 0; }
cat "$LOG"
fails=$(sed -n 's/^ *---- adapter checks, \([0-9][0-9]*\) failed$/\1/p' "$LOG")
rm -f "$LOG"
[ "${fails:-1}" -eq 0 ]
