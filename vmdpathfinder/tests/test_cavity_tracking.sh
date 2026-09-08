#!/bin/sh
# Cavities must keep their identity across frames - see cavity_tracking_check.tcl
# for why a per-frame rank cannot serve as one.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
echo "=============================================================="
echo "cavity-tracking: $DIR/../vmdpathfinder.tcl"
pass=0; fail=0
if tclsh "$DIR/cavity_tracking_check.tcl"; then
    pass=$((pass+1))
else
    fail=$((fail+1))
fi
echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
