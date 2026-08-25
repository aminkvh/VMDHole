#!/bin/sh
# Wrapper for hydro_qco_parity.tcl - plain tclsh, no VMD needed.
here=$(cd "$(dirname "$0")" && pwd)
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
VMDHOLE_TEST_DIR="$here" tclsh "$here/hydro_qco_parity.tcl"
