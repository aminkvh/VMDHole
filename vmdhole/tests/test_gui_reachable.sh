#!/bin/sh
# Wrapper for gui_reachable.tcl. Needs Tk 8.5+ and a real X display, so it SKIPS
# rather than fails where either is absent - the rest of the suite is headless
# and must stay runnable without a screen.
#
# VMD_BIN can name the interpreter. The default prefers vmd2, but the whole file
# passes under VMD 1.9.4a57 (Tk 8.5.6) too - verified with VMD_BIN=vmd.
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
LOG="${TMPDIR:-/tmp}/vmdhole_gui_reachable.$$.log"

# Run on a VIRTUAL display by default. VMD maps real windows, so running this
# on the user's own :0/:1 pops a window up on their desktop every invocation -
# unusable while they are working. Xvfb gives a real X server with no screen.
# VMDHOLE_GUI_DISPLAY=:1 (or any value) forces a specific display; VMDHOLE_GUI_XVFB=0
# opts out entirely and uses whatever DISPLAY is already set.
_XVFB_PID=""
if [ "${VMDHOLE_GUI_XVFB:-1}" != "0" ] && [ -z "${VMDHOLE_GUI_DISPLAY:-}" ] \
        && command -v Xvfb >/dev/null 2>&1; then
    for _d in 99 98 97 96; do
        [ -e "/tmp/.X11-unix/X$_d" ] && continue
        Xvfb ":$_d" -screen 0 1600x1200x24 >/dev/null 2>&1 &
        _XVFB_PID=$!
        # Wait for the socket rather than sleeping a fixed amount.
        for _i in 1 2 3 4 5 6 7 8 9 10; do
            [ -e "/tmp/.X11-unix/X$_d" ] && break
            sleep 0.3
        done
        if [ -e "/tmp/.X11-unix/X$_d" ]; then
            DISPLAY=":$_d"; export DISPLAY
        else
            kill "$_XVFB_PID" 2>/dev/null; _XVFB_PID=""
        fi
        break
    done
fi
[ -n "${VMDHOLE_GUI_DISPLAY:-}" ] && { DISPLAY="$VMDHOLE_GUI_DISPLAY"; export DISPLAY; }
# Kill the private server on exit however we leave.
if [ -n "$_XVFB_PID" ]; then
    trap 'kill "$_XVFB_PID" 2>/dev/null' EXIT INT TERM
fi

if [ -z "$DISPLAY" ]; then
    echo "SKIP: no DISPLAY and no Xvfb - the GUI reachability check needs an X server"
    exit 0
fi
echo "gui-reachability display: $DISPLAY${_XVFB_PID:+ (private Xvfb)}"
BIN="${VMD_BIN:-}"
if [ -z "$BIN" ]; then
    for c in vmd2 vmd; do command -v "$c" >/dev/null 2>&1 && { BIN=$c; break; }; done
fi
[ -n "$BIN" ] || { echo "SKIP: no vmd on PATH"; exit 0; }

echo "gui-reachability: $BIN"
REF="$DIR/fixtures/mole_reference"
run_gui() {   # $1 pdb, $2 het residue to require, $3 selection, $4 start point
    GUI_TEST_PDB="$1" GUI_TEST_HET="$2" GUI_TEST_SEL="$3" GUI_TEST_START="$4" \
    GUI_TEST_ENGINE="${GUI_TEST_ENGINE:-$DIR/../../native/mole_tunnel_engine}" \
    GUI_TEST_LOG="$LOG" VMDHOLE_TCL="$DIR/../vmdhole.tcl" \
        timeout 240 "$BIN" -e "$DIR/gui_reachable.tcl" >/dev/null 2>&1
}
if [ -n "${GUI_TEST_HET+x}" ]; then export GUI_TEST_HET; fi
run_gui "${GUI_TEST_PDB:-$REF/1BL8.pdb}" "" "" ""
if [ ! -s "$LOG" ]; then
    echo "SKIP: the GUI did not open under $BIN (Tk 8.5+ missing?)"
    rm -f "$LOG"; exit 0
fi
cat "$LOG"
# Second pass on a structure that HAS het residues. 1BL8 has none, so the HET
# row in the Lining window renders untested on it - the row is built from the
# engine's H lines. The start point is the origin of the tunnel that MOLE's
# own 1ERI reference reports carrying DA 6 B and DA 7 B; the GUI needs an
# explicit start point, so auto origins are not an option here.
if [ -z "${GUI_TEST_PDB:-}" ] && [ -f "$REF/1ERI.pdb" ]; then
    HETLOG="$LOG.het"
    LOG_SAVE="$LOG"; LOG="$HETLOG"
    run_gui "$REF/1ERI.pdb" "DA" "protein or nucleic" "23.577 28.862 1.087"
    # VMD sometimes dies mid-pass under load, right after the tunnel-run
    # block (the exit-path crash the comment below the verdict already
    # documents, striking earlier) - the pass is fine when it completes:
    # verified 3x directly and via this wrapper. One retry, SAID out loud,
    # only when the first attempt provably died before finishing.
    if ! grep -q 'ALL CHECKS COMPLETE' "$HETLOG" 2>/dev/null; then
        echo "  NOTE  het pass died mid-run (VMD crash under load) - retrying once"
        : > "$HETLOG"
        run_gui "$REF/1ERI.pdb" "DA" "protein or nucleic" "23.577 28.862 1.087"
    fi
    LOG="$LOG_SAVE"
    if [ -s "$HETLOG" ]; then
        grep -E "HET residue|the run produced tunnels" "$HETLOG"
        cat "$HETLOG" >> "$LOG"
    else
        echo "  SKIP  het pass did not open a GUI"
    fi
    rm -f "$HETLOG"
fi
# Verdict from the report lines, not from a trailing summary: `vmd -e` sometimes
# abandons the file after the last command that touched Tk, and VMD segfaults in
# its own exit path on this file, so the summary line is not reliably reached.
# Completion is proven by an explicit marker emitted after the LAST unconditional
# statement, which makes a truncated run a failure rather than a silent pass.
fails=$(grep -c 'FAIL' "$LOG")
# This marker used to be 'show_tunnel_lining survives an empty result set', which
# resolves at line 143 of 5728 - nearly FIRST, not last, despite the comment
# above claiming otherwise. That is how ~5,580 lines could stop executing (the
# A1/A2 block never ran for weeks) while this group still reported complete.
# One marker per pass, and the het pass appends to the same log, so the count is
# the number of passes that ran to completion - at least one, not exactly one.
done_marker=$(grep -c 'ALL CHECKS COMPLETE' "$LOG")
# The results phase either ran or said why; a log missing both is truncated.
res_marker=$(grep -cE 'the run produced tunnels|results phase not run' "$LOG")
# The het pass is the only thing that can produce this line; if 1ERI.pdb is
# present the pass must have reached it rather than dying silently.
if [ -f "$REF/1ERI.pdb" ] && [ -z "${GUI_TEST_PDB:-}" ]; then
    het_marker=$(grep -c 'shows the HET residue' "$LOG")
else
    het_marker=1
fi
rm -f "$LOG"
[ "$fails" -eq 0 ] && [ "$done_marker" -ge 1 ] && [ "$res_marker" -ge 1 ] \
    && [ "$het_marker" -ge 1 ]
