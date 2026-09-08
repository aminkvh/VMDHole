#!/bin/sh
# Wrapper for mean_axis_check.tcl: the Mean Profile 3D tube must lie along the
# CVECT the run was measured with, not a PCA fit. Needs Tk and a real X display
# (build_and_show_mean_surface is GUI-invoked code), so it runs on a VIRTUAL
# display and SKIPS where one cannot be had - same contract as
# test_gui_reachable.sh, whose Xvfb handling this mirrors.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
LOG="${TMPDIR:-/tmp}/vmdpathfinder_mean_axis.$$.log"
VMD_BIN="${VMD_BIN:-vmd2}"
command -v "$VMD_BIN" >/dev/null 2>&1 || { echo "SKIP: no vmd on PATH"; exit 0; }

_XVFB_PID=""
if [ "${VMDPATHFINDER_GUI_XVFB:-1}" != "0" ] && [ -z "${VMDPATHFINDER_GUI_DISPLAY:-}" ] \
        && command -v Xvfb >/dev/null 2>&1; then
    for _d in 95 94 93 92; do
        [ -e "/tmp/.X11-unix/X$_d" ] && continue
        Xvfb ":$_d" -screen 0 1600x1200x24 >/dev/null 2>&1 &
        _XVFB_PID=$!
        for _i in 1 2 3 4 5 6 7 8 9 10; do
            [ -e "/tmp/.X11-unix/X$_d" ] && break
            sleep 0.3
        done
        [ -e "/tmp/.X11-unix/X$_d" ] && { DISPLAY=":$_d"; export DISPLAY; }
        break
    done
fi
[ -n "${VMDPATHFINDER_GUI_DISPLAY:-}" ] && { DISPLAY="$VMDPATHFINDER_GUI_DISPLAY"; export DISPLAY; }
[ -n "${DISPLAY:-}" ] || { echo "SKIP: no DISPLAY and no Xvfb"; exit 0; }

VMDPATHFINDER_TEST_DIR="$DIR"; export VMDPATHFINDER_TEST_DIR
echo "mean-axis: $VMD_BIN on $DISPLAY"
# vmd -e quits the moment stdin reaches EOF, and this script calls [update],
# which is exactly where VMD notices. Standalone that never fires (stdin is the
# terminal); under the suite stdin is already at EOF, so VMD exited right after
# the GUI opened and the group failed ONLY inside the suite. Inheriting stdin is
# therefore not enough - feed it a pipe nothing ever closes, so the outcome does
# not depend on how the caller was invoked.
FIFO="$LOG.fifo"
mkfifo "$FIFO" 2>/dev/null || FIFO=""
if [ -n "$FIFO" ]; then
    sleep 900 > "$FIFO" &
    _SLEEP_PID=$!
    timeout 600 "$VMD_BIN" -e "$DIR/mean_axis_check.tcl" > "$LOG" 2>&1 < "$FIFO"
    rc=$?
    kill "$_SLEEP_PID" 2>/dev/null
    rm -f "$FIFO"
else
    timeout 600 "$VMD_BIN" -e "$DIR/mean_axis_check.tcl" > "$LOG" 2>&1
    rc=$?
fi
[ -n "$_XVFB_PID" ] && kill "$_XVFB_PID" 2>/dev/null
grep -E "  (PASS|FAIL)( |$)|^SKIP:|^  ----" "$LOG"
if grep -q "^SKIP:" "$LOG"; then rm -f "$LOG"; exit 0; fi
if ! grep -q -- "---- mean-axis checks" "$LOG"; then
    echo "FAIL: the check did not run to completion - full log follows"
    tail -25 "$LOG"; exit 1
fi
# Keep the log when a check FAILS - deleting it is what made the first in-suite
# failure undiagnosable.
grep -q -- "---- mean-axis checks, 0 failed" "$LOG" || { echo "  (log kept: $LOG)"; exit 1; }
rm -f "$LOG"
exit 0
