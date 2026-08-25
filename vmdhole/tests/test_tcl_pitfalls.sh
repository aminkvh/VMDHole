#!/bin/sh
# Tcl-syntax pitfalls that parse fine at source time but throw at CALL time.
#
# The defect this exists for: a `switch` body is parsed as a LIST, so a '#' line
# between its cases is NOT a comment - its words become spurious pattern/body
# pairs, and a quoted phrase followed by a comma is invalid list syntax outright.
# property_meta carried such a comment between its gz and dens cases and threw
#   list element in quotes followed by "," instead of space
# on EVERY call, for every property, from the commit that added the comment until
# it was found by a user running the Property heat map. Nothing caught it: the
# file sources cleanly, `info complete` is happy, and no existing test called the
# proc. Sourcing a file proves nothing about its switch bodies.
#
# Comments INSIDE a case's braced body are fine and common - only top-level lines
# between cases are the bug, so that is what this checks.

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SRC="$DIR/../vmdhole.tcl"
echo "=============================================================="
echo "tcl-pitfalls: $SRC"
pass=0; fail=0
ok()   { echo "  PASS  $1"; pass=$((pass+1)); }
bad()  { echo "  FAIL  $1"; fail=$((fail+1)); }

[ -f "$SRC" ] || { echo "  FAIL  source not found"; exit 1; }

# ---- 1. static lint: defects that source cleanly and throw when called ----
if out=$(python3 "$DIR/tcl_lint.py" "$SRC" 2>&1); then
    ok "static lint clean (switch-comment, switch-parity, dead-callback, undeclared-var)"
    echo "$out" | sed 's/^/       /'
else
    bad "static lint found defects:"
    echo "$out" | sed 's/^/       /'
fi

# ---- 2. runtime: every switch-dispatch metadata proc actually returns ----
VMD=${VMD:-vmd}
# The static lint above needs only tclsh; this half needs a real VMD. Without
# one it must SKIP (per-assertion form - the group still ran its static half),
# not FAIL: a bare machine / CI runner has no VMD by design.
if ! command -v "$VMD" >/dev/null 2>&1; then
    echo "  SKIP  runtime metadata-switch check (no vmd on PATH)"
    echo "  -> $pass passed, $fail failed"
    [ "$fail" -eq 0 ]; exit $?
fi
LOG=$(mktemp)
cat > "$LOG.tcl" <<'EOF'
# Source the file UNDER TEST explicitly. VMD auto-loads the DEPLOYED plugin at
# startup, so without this the runtime half silently exercises ~/.vmd/plugins
# instead of the dev tree - it stayed green against a file with the bug
# reintroduced, which is exactly the failure mode this suite exists to prevent.
package provide Tk 8.5
if {[catch {source $::env(PITFALL_SRC)} e]} {
    set out [open "$::env(PITFALL_OUT)" w]
    puts $out "FAIL cannot source file under test: $e"
    puts $out "BAD=1"
    close $out
    quit
}
set out [open "$::env(PITFALL_OUT)" w]
set bad 0
foreach p {kd ww kr charge polarity lipophilicity esp gz dens pfdens zzz-unknown} {
    if {[catch {::VMDHole::property_meta $p} r]} { puts $out "FAIL property_meta $p: $r"; incr bad }
}
foreach k {kd ww kr charge polarity lipophilicity esp gz dens pfdens zzz-unknown} {
    if {[catch {::VMDHole::scheme_display_label $k}]} { puts $out "FAIL scheme_display_label $k"; incr bad }
}
puts $out "BAD=$bad"
close $out
quit
EOF
PITFALL_OUT="$LOG.out" PITFALL_SRC="$SRC" $VMD -dispdev text -e "$LOG.tcl" >/dev/null 2>&1
if [ -f "$LOG.out" ] && grep -q "^BAD=0$" "$LOG.out"; then
    ok "property_meta / scheme_display_label return for every scheme (+unknown)"
else
    bad "a metadata switch threw when called:"
    [ -f "$LOG.out" ] && sed 's/^/          /' "$LOG.out"
fi
rm -f "$LOG" "$LOG.tcl" "$LOG.out"

echo "  -> $pass passed, $fail failed"
[ "$fail" -eq 0 ]
