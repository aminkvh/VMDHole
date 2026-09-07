#!/bin/sh
# Wrapper for hole_tcl_fallback.tcl - plain tclsh, no VMD needed.
here=$(cd "$(dirname "$0")" && pwd)
command -v tclsh >/dev/null 2>&1 || { echo "SKIP: no tclsh"; exit 0; }
tclsh "$here/hole_tcl_fallback.tcl"
